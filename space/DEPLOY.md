# Deploying COAL to a free Hugging Face Space (not done yet)

Everything is prepared; nobody has deployed it. It takes about 5 minutes.

1. **Check the Gradio pin first.** `space/README.md` (`sdk_version`), `space/requirements.txt` and the
   `gui` extra in `pyproject.toml` all say `5.49.1`. That pin was written without PyPI access. Run the GUI
   locally (`pip install -e .[gui]`, `python -m coal_app`); if you change the version, change all three.
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
7. **MCP for AI agents**: the Space serves MCP at
   `https://<you>-coal.hf.space/gradio_api/mcp/sse` (the Space's "Use via API / MCP" link shows the exact
   config).

Notes
- Free CPU Spaces sleep after inactivity; the first visit then takes ~30 s to wake up.
- Uploaded files and results live in the Space's temporary folder and vanish on restart.
- Everything runs on the Space's CPU; large inputs (hundreds of MB) will be slow there. Point heavy users
  to `pip install` and the command line.
- To update later: re-run step 4's `build_space.py` into the same clone, then commit and push.
