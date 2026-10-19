"""Verified native-binary reuse must reject changed or failed build evidence."""
import importlib.util
import json
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location(
    'local_pic_runner', Path(__file__).parents[1]/'research/pic-local/run_local_pulsar.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


@pytest.fixture
def built(tmp_path):
    binary, pgen, receipt = [tmp_path/name for name in ('native.x', 'pgen.hpp', 'verification.json')]
    binary.write_bytes(b'verified executable bytes')
    pgen.write_bytes(b'verified external source')
    sources = [{'repository': 'pinned-upstream', 'commit': 'a'*40}]
    content = {'exit_code': 0, 'build_commands': [['configure'], ['compile']],
               'source_pins': sources, 'external_pgen_sha256': runner.digest(pgen),
               'native_executable_sha256': runner.digest(binary)}
    receipt.write_text(json.dumps(content))
    return binary, pgen, receipt, sources, content


def test_exact_completed_build_can_be_reused(built):
    binary, pgen, receipt, sources, _ = built
    record = runner.verify_build_receipt(receipt, tuple(sources), binary, pgen)
    assert record['sha256'] == runner.digest(receipt)


def test_changed_native_binary_is_rejected(built):
    binary, pgen, receipt, sources, _ = built
    binary.write_bytes(b'changed executable')
    with pytest.raises(ValueError, match='Native binary differs'):
        runner.verify_build_receipt(receipt, sources, binary, pgen)


def test_changed_problem_generator_is_rejected(built):
    binary, pgen, receipt, sources, _ = built
    pgen.write_bytes(b'changed external source')
    with pytest.raises(ValueError, match='Problem generator differs'):
        runner.verify_build_receipt(receipt, sources, binary, pgen)


def test_failed_or_chained_receipts_and_different_pins_are_rejected(built):
    binary, pgen, receipt, sources, content = built
    for changes in ({'exit_code': 1}, {'build_commands': []},
                    {'reused_verified_build': {'path': 'another.json'}},
                    {'source_pins': [{'commit': 'b'*40}]}):
        receipt.write_text(json.dumps({**content, **changes}))
        with pytest.raises(ValueError, match='direct successful build'):
            runner.verify_build_receipt(receipt, sources, binary, pgen)
