#ifndef PROBLEM_GENERATOR_H
#define PROBLEM_GENERATOR_H

// CPU-only local kinetic prototype. This external problem generator leaves the
// pinned Entity source unchanged. There is no atmosphere/radiation prediction.
#include "enums.h"
#include "global.h"
#include "arch/kokkos_aliases.h"
#include "traits/pgen.h"
#include "framework/domain/metadomain.h"
#include "framework/parameters/parameters.h"
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace user {
using namespace ntt;

template<Dimension D> struct Dipole {
  real_t radius, mu, chi, h;
  Dipole(real_t r, real_t bp, real_t c, real_t dx)
    : radius(r), mu(bp*r*r*r/2), chi(c), h(dx) {}
  Inline real_t potential(int c, const coord_t<D>& x) const {
    const real_t r2=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
    const real_t f=r2>=radius*radius ? 1/(math::sqrt(r2)*r2)
      : (2.5-1.5*r2/(radius*radius))/(radius*radius*radius);
    const real_t mx=mu*math::sin(chi), mz=mu*math::cos(chi);
    return (c==0 ? -mz*x[1] : c==1 ? mz*x[0]-mx*x[2] : mx*x[1])*f;
  }
  Inline real_t derivative(int c, int d, const coord_t<D>& x) const {
    coord_t<D> xp{ZERO},xm{ZERO};
    for(unsigned int j=0;j<D;++j){xp[j]=x[j];xm[j]=x[j];}
    xp[d]+=h/2; xm[d]-=h/2;
    return (potential(c,xp)-potential(c,xm))/h;
  }
  Inline real_t bx1(const coord_t<D>& x) const {return derivative(2,1,x)-derivative(1,2,x);}
  Inline real_t bx2(const coord_t<D>& x) const {return derivative(0,2,x)-derivative(2,0,x);}
  Inline real_t bx3(const coord_t<D>& x) const {return derivative(1,0,x)-derivative(0,1,x);}
  Inline real_t ex1(const coord_t<D>&) const {return 0;}
  Inline real_t ex2(const coord_t<D>&) const {return 0;}
  Inline real_t ex3(const coord_t<D>&) const {return 0;}
};

template<Dimension D> struct ZeroFields {
  Inline real_t bx1(const coord_t<D>&) const {return 0;}
  Inline real_t bx2(const coord_t<D>&) const {return 0;}
  Inline real_t bx3(const coord_t<D>&) const {return 0;}
  Inline real_t ex1(const coord_t<D>&) const {return 0;}
  Inline real_t ex2(const coord_t<D>&) const {return 0;}
  Inline real_t ex3(const coord_t<D>&) const {return 0;}
};

template<SimEngine::type S,class M> struct PGen {
  static constexpr auto D=M::Dim;
  static constexpr auto engines=::traits::pgen::compatible_with<SimEngine::SRPIC>{};
  static constexpr auto metrics=::traits::pgen::compatible_with<Metric::Minkowski>{};
  static constexpr auto dimensions=::traits::pgen::compatible_with<Dim::_3D>{};
  const SimulationParams& params;
  const real_t R, Bpolar, chi, omega, spinup, h, dt, B0, V0, q0, n0, vout;
  const int pairs, stride;
  const unsigned int maxsteps;
  const std::string directory;
  Dipole<D> init_flds;
  std::ofstream events, ledger;
  std::vector<std::set<npart_t>> recorded_dead;
  double conductor_work=0, injected_energy=0, star_lost_energy=0, boundary_lost_energy=0;
  unsigned long long injected_pairs=0;

  PGen(const SimulationParams& p,const Metadomain<S,M>& m)
    : params(p), R(p.template get<real_t>("setup.star_radius",.5)),
      Bpolar(p.template get<real_t>("setup.polar_field",1)),
      chi(p.template get<real_t>("setup.inclination_deg",75)*M_PI/180),
      omega(p.template get<real_t>("setup.omega",1)),
      spinup(p.template get<real_t>("setup.spinup_time",.8)),
      h(m.mesh().metric.get_dx()),
      dt(p.template get<real_t>("algorithms.timestep.CFL")*m.mesh().metric.find_dxMin()),
      B0(p.template get<real_t>("scales.B0")), V0(p.template get<real_t>("scales.V0")),
      q0(p.template get<real_t>("scales.q0")), n0(p.template get<real_t>("scales.n0")),
      vout(p.template get<real_t>("setup.radial_speed",.2)),
      pairs(p.template get<int>("setup.pairs_per_step",64)),
      stride(p.template get<int>("setup.snapshot_interval",20)),
      maxsteps(static_cast<unsigned int>(p.template get<real_t>("simulation.runtime")/dt)),
      directory(p.template get<std::string>("setup.output_directory","native")),
      init_flds(R,Bpolar,chi,h), recorded_dead(2) {
    if(R<=2*h || R*omega>=.8 || pairs<=0 || stride<=0 || spinup<=0)
      throw std::runtime_error("Invalid resolved source/injection configuration");
    std::filesystem::create_directories(directory);
    events.open(directory+"/events.csv",std::ios::out|std::ios::trunc);
    ledger.open(directory+"/ledger.csv",std::ios::out|std::ios::trunc);
    events<<"completed_step,x_time,u_time,event,species,id,x,y,z,ux,uy,uz,weight\n";
    ledger<<"completed_step,time,field_energy,particle_total_energy,particle_kinetic_energy,"
            "conductor_work_cumulative,injected_energy_cumulative,star_loss_energy_cumulative,"
            "boundary_loss_energy_cumulative,alive_electrons,alive_positrons,divB_rms_active,"
            "gauss_rms_active,injected_pairs\n";
    events<<std::setprecision(17); ledger<<std::setprecision(17);
  }

  auto MatchFields(simtime_t) const -> ZeroFields<D> {return {};}
  real_t angular_velocity(real_t t) const {
    const real_t s=std::min(real_t(1),std::max(real_t(0),t/spinup));
    return omega*s*s*(3-2*s);
  }
  real_t spin_phase(real_t t) const {
    const real_t s=std::min(real_t(1),std::max(real_t(0),t/spinup));
    return omega*(spinup*(s*s*s-.5*s*s*s*s)+std::max(real_t(0),t-spinup));
  }
  void event(unsigned int step,real_t t,const char* kind,unsigned int sp,
             npart_t id,real_t x,real_t y,real_t z,real_t ux,real_t uy,real_t uz,real_t w) {
    events<<step<<','<<t<<','<<t-dt/2<<','<<kind<<','<<sp<<','<<id<<','
      <<x<<','<<y<<','<<z<<','<<ux<<','<<uy<<','<<uz<<','<<w<<'\n';
  }

  void inject(Domain<S,M>& dom,unsigned int step,real_t t) {
    const real_t radius=R+.7*h, om=angular_velocity(t);
    for(unsigned int s=0;s<2;++s) {
      auto& sp=dom.species[s]; const auto begin=sp.npart(), counter=sp.counter();
      if(begin+pairs>=sp.maxnpart()) throw std::runtime_error("Particle capacity exhausted");
      for(unsigned int j=0;j<pairs;++j) {
        const real_t zunit=1-2*(real_t(j)+.5)/pairs;
        const real_t az=2.399963229728653*(j+real_t(step)*pairs);
        const real_t rxy=std::sqrt(1-zunit*zunit);
        const real_t x=radius*rxy*std::cos(az), y=radius*rxy*std::sin(az), z=radius*zunit;
        const real_t vx=vout*x/radius-om*y, vy=vout*y/radius+om*x, vz=vout*z/radius;
        const real_t gamma=1/std::sqrt(1-vx*vx-vy*vy-vz*vz);
        const real_t ux=gamma*vx, uy=gamma*vy, uz=gamma*vz;
        coord_t<D> cd{ZERO}; dom.mesh.metric.template convert<Crd::Ph,Crd::Cd>({x,y,z},cd);
        const auto n=begin+j;
        sp.i1(n)=static_cast<int>(cd[0]); sp.dx1(n)=cd[0]-sp.i1(n);
        sp.i2(n)=static_cast<int>(cd[1]); sp.dx2(n)=cd[1]-sp.i2(n);
        sp.i3(n)=static_cast<int>(cd[2]); sp.dx3(n)=cd[2]-sp.i3(n);
        sp.i1_prev(n)=sp.i1(n); sp.i2_prev(n)=sp.i2(n); sp.i3_prev(n)=sp.i3(n);
        sp.dx1_prev(n)=sp.dx1(n); sp.dx2_prev(n)=sp.dx2(n); sp.dx3_prev(n)=sp.dx3(n);
        sp.ux1(n)=ux; sp.ux2(n)=uy; sp.ux3(n)=uz; sp.weight(n)=1;
        sp.tag(n)=ParticleTag::alive; sp.pld_i(n,pldi::spcCtr)=counter+j;
        injected_energy+=q0*gamma;
        coord_t<D> stored{ZERO};
        dom.mesh.metric.template convert<Crd::Cd,Crd::Ph>(
          {static_cast<real_t>(sp.i1(n))+static_cast<real_t>(sp.dx1(n)),
           static_cast<real_t>(sp.i2(n))+static_cast<real_t>(sp.dx2(n)),
           static_cast<real_t>(sp.i3(n))+static_cast<real_t>(sp.dx3(n))},stored);
        event(step,t,"birth",s+1,counter+j,stored[0],stored[1],stored[2],ux,uy,uz,1);
      }
      sp.set_npart(begin+pairs); sp.set_counter(counter+pairs);
    }
    injected_pairs+=pairs;
  }

  void absorb(Domain<S,M>& dom,unsigned int step,real_t t) {
    for(unsigned int s=0;s<2;++s) {
      auto& sp=dom.species[s];
      for(npart_t p=0;p<sp.npart();++p) {
        const auto id=sp.pld_i(p,pldi::spcCtr);
        if(recorded_dead[s].contains(id)) continue;
        coord_t<D> x{ZERO};
        dom.mesh.metric.template convert<Crd::Cd,Crd::Ph>(
          {static_cast<real_t>(sp.i1(p))+static_cast<real_t>(sp.dx1(p)),
           static_cast<real_t>(sp.i2(p))+static_cast<real_t>(sp.dx2(p)),
           static_cast<real_t>(sp.i3(p))+static_cast<real_t>(sp.dx3(p))},x);
        const bool inside=x[0]*x[0]+x[1]*x[1]+x[2]*x[2]<R*R;
        const bool dead=sp.tag(p)!=ParticleTag::alive;
        if(inside || dead) {
          const real_t u2=sp.ux1(p)*sp.ux1(p)+sp.ux2(p)*sp.ux2(p)+sp.ux3(p)*sp.ux3(p);
          const double loss=q0*sp.weight(p)*std::sqrt(1+u2);
          if(dead) boundary_lost_energy+=loss; else star_lost_energy+=loss;
          event(step,t,dead?"boundary_death":"star_death",s+1,id,x[0],x[1],x[2],
                sp.ux1(p),sp.ux2(p),sp.ux3(p),sp.weight(p));
          recorded_dead[s].insert(id); sp.tag(p)=ParticleTag::dead;
        }
      }
    }
  }

  void conductor(Domain<S,M>& dom,real_t t) {
    auto F=dom.fields.em; const auto N=dom.mesh.n_active(in::x1); const int g=N_GHOSTS;
    const real_t om=angular_velocity(t), x0=dom.mesh.extent(in::x1).first;
    double work=0;
    for(int i=g;i<g+int(N);++i) for(int j=g;j<g+int(N);++j)
      for(int k=g;k<g+int(N);++k) {
        const real_t x=x0+(i-g)*h,y=x0+(j-g)*h,z=x0+(k-g)*h;
        if((x+h/2)*(x+h/2)+y*y+z*z<R*R) {
          const real_t target=-om*(x+h/2)*.5*(F(i,j,k,em::bx3)+F(i,j-1,k,em::bx3));
          work+=target*target-F(i,j,k,em::ex1)*F(i,j,k,em::ex1); F(i,j,k,em::ex1)=target;
        }
        if(x*x+(y+h/2)*(y+h/2)+z*z<R*R) {
          const real_t target=-om*(y+h/2)*.5*(F(i,j,k,em::bx3)+F(i-1,j,k,em::bx3));
          work+=target*target-F(i,j,k,em::ex2)*F(i,j,k,em::ex2); F(i,j,k,em::ex2)=target;
        }
        if(x*x+y*y+(z+h/2)*(z+h/2)<R*R) {
          const real_t target=om*(y*.5*(F(i,j,k,em::bx2)+F(i-1,j,k,em::bx2))+
                                      x*.5*(F(i,j,k,em::bx1)+F(i,j-1,k,em::bx1)));
          work+=target*target-F(i,j,k,em::ex3)*F(i,j,k,em::ex3); F(i,j,k,em::ex3)=target;
        }
      }
    conductor_work+=.5*B0*B0*h*h*V0*work;
  }

  void dump(Domain<S,M>& dom,unsigned int step,real_t t) {
    Kokkos::fence();
    const int n=dom.mesh.n_active(in::x1),g=N_GHOSTS,a=n+2*g;
    const auto F=dom.fields.em;
    const auto J=dom.fields.cur;
    std::vector<double> rho(std::size_t(a)*a*a,0);
    auto index=[a](int i,int j,int k){return (std::size_t(i)*a+j)*a+k;};
    std::ostringstream ss; ss<<std::setw(6)<<std::setfill('0')<<step;
    const std::string stem="state-"+ss.str();
    std::ofstream particles(directory+"/"+stem+"-particles.csv");
    particles<<std::setprecision(17)<<"species,id,x,y,z,ux,uy,uz,weight\n";
    unsigned long long alive[2]={0,0}; double kinetic=0,total=0;
    for(unsigned int s=0;s<2;++s) {
      const auto& sp=dom.species[s];
      for(npart_t p=0;p<sp.npart();++p) if(sp.tag(p)==ParticleTag::alive) {
        ++alive[s]; coord_t<D> x{ZERO};
        dom.mesh.metric.template convert<Crd::Cd,Crd::Ph>(
          {static_cast<real_t>(sp.i1(p))+static_cast<real_t>(sp.dx1(p)),
           static_cast<real_t>(sp.i2(p))+static_cast<real_t>(sp.dx2(p)),
           static_cast<real_t>(sp.i3(p))+static_cast<real_t>(sp.dx3(p))},x);
        const double gamma=std::sqrt(1+sp.ux1(p)*sp.ux1(p)+sp.ux2(p)*sp.ux2(p)+sp.ux3(p)*sp.ux3(p));
        total+=q0*gamma*sp.weight(p); kinetic+=q0*(gamma-1)*sp.weight(p);
        particles<<s+1<<','<<sp.pld_i(p,pldi::spcCtr)<<','<<x[0]<<','<<x[1]<<','<<x[2]<<','
          <<sp.ux1(p)<<','<<sp.ux2(p)<<','<<sp.ux3(p)<<','<<sp.weight(p)<<'\n';
        const int ii=sp.i1(p)+g,jj=sp.i2(p)+g,kk=sp.i3(p)+g;
        const double f[3]={sp.dx1(p),sp.dx2(p),sp.dx3(p)};
        for(int u=0;u<2;++u) for(int v=0;v<2;++v) for(int w=0;w<2;++w) {
          const int i=ii+u,j=jj+v,k=kk+w;
          if(i>=0&&i<a&&j>=0&&j<a&&k>=0&&k<a)
            rho[index(i,j,k)]+=sp.charge()*sp.weight(p)*(u?f[0]:1-f[0])*
              (v?f[1]:1-f[1])*(w?f[2]:1-f[2])*q0/(B0*V0);
        }
      }
    }
    std::vector<double> out(std::size_t(a)*a*a*10);
    double field_energy=0,divb2=0,gauss2=0; std::size_t measured=0;
    const real_t origin=dom.mesh.extent(in::x1).first;
    for(int i=0;i<a;++i) for(int j=0;j<a;++j) for(int k=0;k<a;++k) {
      const auto ix=index(i,j,k); for(int c=0;c<6;++c)out[ix*10+c]=F(i,j,k,c);
      for(int c=0;c<3;++c)out[ix*10+6+c]=J(i,j,k,c); out[ix*10+9]=rho[ix];
      for(int c=0;c<10;++c)if(!std::isfinite(out[ix*10+c]))throw std::runtime_error("Nonfinite exported field");
      if(i>=g&&i<g+n&&j>=g&&j<g+n&&k>=g&&k<g+n) {
        for(int c=0;c<6;++c)field_energy+=.5*B0*B0*h*h*V0*F(i,j,k,c)*F(i,j,k,c);
        const double x=origin+(i-g)*h,y=origin+(j-g)*h,z=origin+(k-g)*h;
        const double r=std::sqrt(x*x+y*y+z*z);
        if(r>R+3*h&&r<1.8) {
          const double divb=F(i+1,j,k,em::bx1)-F(i,j,k,em::bx1)+
            F(i,j+1,k,em::bx2)-F(i,j,k,em::bx2)+F(i,j,k+1,em::bx3)-F(i,j,k,em::bx3);
          const double dive=F(i,j,k,em::ex1)-F(i-1,j,k,em::ex1)+
            F(i,j,k,em::ex2)-F(i,j-1,k,em::ex2)+F(i,j,k,em::ex3)-F(i,j,k-1,em::ex3);
          divb2+=divb*divb; gauss2+=(dive-rho[ix])*(dive-rho[ix]); ++measured;
        }
      }
    }
    std::ofstream binary(directory+"/"+stem+"-fields.bin",std::ios::binary);
    binary.write(reinterpret_cast<const char*>(out.data()),out.size()*sizeof(double));
    std::ofstream meta(directory+"/"+stem+".json"); meta<<std::setprecision(17);
    meta<<"{\n\"completed_step\":"<<step<<",\n\"position_and_E_time\":"<<t
      <<",\n\"momentum_B_and_J_time\":"<<t-dt/2<<",\n\"dt\":"<<dt
      <<",\n\"field_shape\":["<<a<<','<<a<<','<<a<<",10],\n\"ghost_cells_each_side\":"<<g
      <<",\n\"active_resolution\":["<<n<<','<<n<<','<<n<<"],\n\"origin\":["<<origin<<','<<origin<<','<<origin
      <<"],\n\"dx\":"<<h<<",\n\"B0\":"<<B0<<",\n\"q0\":"<<q0<<",\n\"n0\":"<<n0<<",\n\"V0\":"<<V0
      <<",\n\"dtype\":\"little-endian float64\",\n\"storage_order\":\"C: i,j,k,component\",\n"
        "\"field_components\":[\"Ex_U\",\"Ey_U\",\"Ez_U\",\"Bx_U\",\"By_U\",\"Bz_U\",\"Jx_U_over_ppc\",\"Jy_U_over_ppc\",\"Jz_U_over_ppc\",\"rho_vertex_gauss_units\"],\n"
        "\"fields_file\":\""<<stem<<"-fields.bin\",\n\"particles_file\":\""<<stem<<"-particles.csv\",\n"
        "\"physical_tetrad_normalized_fields\":\"multiply stored E_U and B_U by dx; multiply again by B0 for code EM units\",\n"
        "\"rho_definition\":\"CIC vertex sum(q_sign*weight)*q0/(B0*V0), compared with backward differences of raw E_U\",\n"
        "\"units\":\"scaled dimensionless kinetic model; c=Omega_final=1; no physical Vela normalization\",\n"
        "\"particle_momentum\":\"spatial proper velocity u=gamma*v/c at staggered half step\",\n"
        "\"stellar_radius\":"<<R<<",\n\"inclination_deg\":"<<chi*180/M_PI
      <<",\n\"spinup_time\":"<<spinup<<",\n\"final_omega\":"<<omega
      <<",\n\"instantaneous_omega\":"<<angular_velocity(t)<<",\n\"integrated_spin_phase\":"<<spin_phase(t)
      <<",\n\"larmor0\":"<<params.template get<real_t>("scales.larmor0")
      <<",\n\"skindepth0\":"<<params.template get<real_t>("scales.skindepth0")<<",\n"
        "\"newborn_momentum_initialization\":\"prescribed half-step injection momenta; no birth-time half-kick or thermal atmosphere solve\",\n"
        "\"radiation_or_optical_prediction\":false,\n\"steady_state_or_convergence_claim\":false\n}\n";
    ledger<<step<<','<<t<<','<<field_energy<<','<<total<<','<<kinetic<<','<<conductor_work<<','
      <<injected_energy<<','<<star_lost_energy<<','<<boundary_lost_energy<<','<<alive[0]<<','<<alive[1]<<','
      <<std::sqrt(divb2/std::max(std::size_t(1),measured))<<','
      <<std::sqrt(gauss2/std::max(std::size_t(1),measured))<<','<<injected_pairs<<'\n';
    ledger.flush(); events.flush();
  }

  void InitPrtls(Domain<S,M>& dom) { inject(dom,0,0); dump(dom,0,0); }
  void CustomPostStep(timestep_t step,simtime_t time,Domain<S,M>& dom) {
    Kokkos::fence(); const auto completed=static_cast<unsigned int>(step)+1;
    const real_t t=time+dt;
    absorb(dom,completed,t); conductor(dom,t); inject(dom,completed,t);
    if(completed%stride==0 || completed>=maxsteps) dump(dom,completed,t);
  }
};
} // namespace user
#endif
