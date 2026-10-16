# v1 package staging notes

The source project contains large local runtimes, build products and historical experiments. Publish a deliberate release selection rather than uploading the entire working directory.

- Include the final 32-second Reel, poster, `instagram-description.text`, `MODEL-NOTES.md`, manifest, and the documented reproduction entry points
- Include the source modules used by the final scene, their scientific tests, native build scripts, configuration, complete upstream patch chain and GPL-2.0 notices
- Preserve the upstream commit identifier and source/binary/snapshot hashes; do not present the modified native solver as unmodified upstream COFFEE
- Include compact final diagnostic reports and local input-data provenance. Keep historical reports clearly named, or describe them separately from the final selection
- Avoid `.venv`, dependency/runtime caches, local native binaries, compiled products, temporary previews and unrelated local experiments
- Large final field snapshots or source data may be supplied as release assets or reproducibly fetched inputs, with exact checksums. All source needed to rebuild any distributed GPL native binary must remain available under its applicable terms
- The caption links the selected repository, `https://github.com/Loganpis/sideratlas-pulsar`, and contains 1,970 characters, including line breaks. Recheck the count after any edits and update the final manifest after all media/notes are frozen

This file records packaging guidance; it is not a statement that the repository or release has already been published.
