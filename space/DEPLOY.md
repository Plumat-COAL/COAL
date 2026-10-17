# Deploying COAL to a free Hugging Face Space (not done yet)

Everything is prepared; nobody has deployed it. It takes about 5 minutes.

1. **The Gradio pin.** `space/README.md` (`sdk_version`), `space/requirements.txt` and the `gui` extra in
   `pyproject.toml` all say `5.49.1`. CI's `gui` job installs `gradio[mcp]==5.49.1` and passes `pip check` and
   the GUI tests (verified by CI on Linux). Still run the GUI locally once (`pip install -e .[gui]`,
   `python -m coal_app`); if you change the version, change all three.
   <!-- CHANGED(R2-D10): was "That pin was written without PyPI access". See overnight/integration/CHANGES.md -->
2. **Create the Space**: https://huggingface.co/new-space -> owner + name `coal`, licence MIT,
   SDK **Gradio**, template **Blank**, hardware **CPU basic (free)**, visibility Public. Create.
3. **Get a write token**: https://huggingface.co/settings/tokens -> New token, role *Write*.
4. **Clone the empty Space and fill it** (from the COAL repo root):
   ```
   git clone https://huggingface.co/spaces/<you>/coal ../coal-space
   python space/build_space.py ../coal-space
   ```
5. **Push** (git asks for a username and password: use your HF username and the token as password):
   ```
   cd ../coal-space
   git add -A
   git commit -m "COAL Space"
   git push
   ```
6. **Watch the build** on the Space page (Logs tab). When it says *Running*, open the App tab and press
   **Try an example -> One broken gene**.
7. **MCP is off on the Space on purpose.** The tools take file paths on the machine that runs them; a remote
   agent cannot give the Space its own files, and on a public server the tools must not read arbitrary
   paths. `space_app.py` sets `COAL_HOSTED=1`, which (a) shows a "your files are processed on this server"
   notice instead of "nothing leaves your computer", (b) restricts `codon_align` / `explain_report` (also
   reachable through Gradio's HTTP API) to the bundled examples, the Space's own run folder and Gradio's
   upload folder (symlinks and `..` are resolved first) and (c) launches without the MCP server. Agents
   should run COAL locally (`python -m coal_app --mcp`).

Notes
- Free CPU Spaces sleep after inactivity; the first visit then takes ~30 s to wake up.
- Uploaded files and results live in the Space's temporary folder; run folders older than 24 hours are
  deleted (at startup and before each new run), and everything vanishes on restart.
- Everything runs on the Space's CPU; large inputs (hundreds of MB) will be slow there. Point heavy users
  to `pip install` and the command line.
- To update later: re-run step 4's `build_space.py` into the same clone, then commit and push.
