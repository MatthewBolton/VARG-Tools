# User guide source

The modular `.qmd` files in this directory are the editable source for the VARG-Tools guide. Do not edit `docs/VARG-Tools_User_Guide.md` or `docs/guide/VARG-Tools_User_Guide.qmd` directly. They are generated for the in-app/standalone HTML guide and the combined Quarto/PDF guide, respectively.

From the app repository root:

```powershell
Rscript docs/guide/build_workflow_map.R
Rscript docs/build_in_app_guide.R
Rscript docs/validate_guide_render.R
```

Run the first command only when the workflow-map source changes. `docs/sync_version.R` is retained as a compatibility wrapper and now rebuilds the full in-app guide.

When Quarto is installed, render the website from this directory with `quarto render`. The in-app guide build does not require Quarto.

Current figure policy:

- prefer concise concept or decision figures over full-window screenshots;
- use each figure once and provide meaningful alt text;
- keep interface labels synchronized with the current source;
- place advanced settings in `reference.qmd`, not in the first-run path.
