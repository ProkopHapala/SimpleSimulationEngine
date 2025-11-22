### Condensed log of our 2-day WebGPU-on-NVIDIA-Linux nightmare (Nov 22 2025)

| Step | What we tried | Result | Insight |
|------|---------------|--------|---------|
| Initial state | Normal Chrome on hybrid laptop | WebGPU = Disabled (blocked on Intel/Mesa) | Chrome blocklists WebGPU on Intel Linux |
| `--enable-features=Vulkan` alone | Worked partially | WebGPU = Hardware accelerated, but still ran on Intel iGPU → demos crashed with “out of device memory” | Vulkan backend enabled, but GPU selection ignored |
| All offload attempts on X11 “on-demand” | `__NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia …` + every possible extra var | Always crashed with “Invalid visual ID requested” EGL error | Long-standing (2019–2025) **unfixed bug** in Chromium + proprietary NVIDIA PRIME offload on X11 |
| `prime-run` | `prime-run google-chrome-stable` | Command not found | `prime-run` missing on many Ubuntu 22.04/24.04 installs (package broken or incomplete) |
| `sudo prime-select nvidia` + reboot | Full switch to NVIDIA | **Worked perfectly** – WebGPU on GTX 1650, boids 200+ fps, full shared memory | Only 100% reliable way on X11 |
| `sudo prime-select intel` + reboot | Back to Intel as default | Cool & long battery again | System-wide switch is safe and reversible |
| Offload after switching back to Intel | `env __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia google-chrome-stable` | **Now works perfectly** | The EGL crash only happens in “on-demand” mode; after explicit `prime-select intel` the offload path becomes clean |

### Final working setup (what you use forever)

1. Run once:  
   `sudo prime-select intel && sudo reboot` → daily life on Intel (cool, long battery)
2. One desktop icon for your molecule editor:  
   ```bash
   env __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia google-chrome-stable
   ```
   → automatically uses GTX 1650 with full WebGPU + real `var<workgroup>` shared memory

Everything else (DRI_PRIME, __VK_LAYER_NV_optimus, --use-gl=desktop, prime-run, etc.) = red herrings in 2025.

You now have the definitive recipe – no one should ever suffer this again.  
Ready for the tiled Gauss–Seidel WGSL code whenever you are! 🚀