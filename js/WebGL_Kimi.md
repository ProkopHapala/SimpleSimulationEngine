
## User 1

I was thinking to program some physical simulation on web using WebGL using user (client-site) GPU. But it seems the support is quite limited, just OpenGL ES, and in GLS shaderts I have problem to use 3D textures.
Also I'm not sure if it is possible to use one shader to pre-calculated some boffers and then use them for rendering by other shader?


---

You're right that WebGL (especially WebGL 1) is quite limited compared to desktop OpenGL, and your two pain points—**3D textures** and **inter-shader data exchange**—are real and well-known.

### ✅ 1. **3D Textures in WebGL**
- **WebGL 1**: No 3D textures at all.
- **WebGL 2**: **Supports 3D textures**, but **you cannot render to a full 3D texture in one pass**. You can only attach **one slice (layer)** of a 3D texture to a framebuffer using `framebufferTextureLayer` .

> This means if you want to simulate a 3D field (like fluid or temperature), you have to **render slice-by-slice**, which is slow and awkward.

### ✅ 2. **Using One Shader to Precompute Buffers for Another**
WebGL **does not allow shaders to write to arbitrary buffers** directly. But you can **simulate this** using **render-to-texture**:

- You **render to a texture** (or multiple textures) using a **framebuffer object (FBO)**.
- Then you **bind those textures** as inputs to the **next shader pass**.

This is the **standard multi-pass rendering pattern** in WebGL. You can chain as many passes as you want, each writing to a texture, and then use those textures in later passes .

> This is how most GPU-based simulations (like fluid solvers or particle systems) are done in WebGL.

### ✅ Summary of Workarounds

| Limitation | Workaround |
|------------|------------|
| No 3D texture writes | Use 2D texture atlas or 2D array textures; render slice-by-slice  |
| No shader-to-buffer writes | Use **render-to-texture** with FBOs; chain multiple shader passes  |

### ✅ Bonus: If You Can, Consider **WebGPU**
If you're targeting modern browsers and want **real 3D textures**, **compute shaders**, and **storage buffers**, **WebGPU** is the future. It supports:
- 3D textures with full read/write
- Compute shaders for simulation
- Proper inter-shader data sharing via storage textures and buffers 

Let me know if you want a minimal WebGL2 example of render-to-texture or a 3D texture simulation loop.

---

## User 2

Ok, ale compute shaders na webgl nejdou nebo jo? Co když chci pomoci GPU simulovat třeba particle systém pomoci molekulární dynamiky (verlet, leapfrog) a modifikovat Vertex-buffery?

---

Compute shadery v **standardním WebGL nefungují** – ani ve WebGL 2.0. Existovala experimentální extenze **WebGL 2.0 Compute**, která je ale dnes **zastaralá a odstraněná z Chromu** (byla dostupná jen pod vlajkou a jen na některých systémech) . Pro GPU výpočty na webu je **WebGPU** správná cesta, ale pokud potřebujete širokou kompatibilitu, musíte použít triky v čistém WebGL 2.

---

### **Simulace částic (Verlet, Leapfrog) ve WebGL 2**

Pro molekulární dynamiku a modifikaci vertex bufferů existují dvě hlavní metody:

#### **1. Transform Feedback (WebGL 2) – NEJLEPŠÍ VOLBA**
Umožňuje **vertex shaderu zapisovat do bufferu**, který pak můžete použít jako vstup pro další kreslicí průchod. To je ideální pro částice:

- **Průchod 1**: Vertex shader počítá nové pozice částic (integrace Verlet/leapfrog) a výstup uloží do bufferu přes `transformFeedback`
- **Průchod 2**: Tentýž buffer (nebo jiný) použijete jako vertex buffer pro vykreslení

**Výhody**: 
- Data zůstávají na GPU, žádný přenos CPU ←→ GPU
- Funguje i pro desítky tisíc částic 

**Limitace**:
- Nejde číst a zároveň zapisovat do stejného bufferu v jednom průchodu (potřebujete ping-pong mezi dvěma buffery)
- Pouze WebGL 2, nefunguje na WebGL 1 

#### **2. Simulace ve fragment shaderu (GPGPU s texturami)**
Pokud potřebujete složitější výpočty:
- Zakódujte stav částic do textury (např. RGBA = position.xyz + velocity.w)
- Použijte fragment shader pro výpočet nového stavu (render-to-texture)
- V dalším průchodu přečtěte texturu ve vertex shaderu (přes `texelFetch`)

**Nevýhoda**: Nutnost manuálního kódování/dekódování dat, horší přesnost (float 32-bit).

---

### **Modifikace vertex bufferů – ano, ale s omezeními**

- **Přes Transform Feedback**: Ano, vertex shader může **zcela přepsat** obsah bufferu 
- **Přes CPU**: Každý frame poslat `bufferSubData()` – pomalé pro velké datasety
- **Přes compute shader**: **NEJDE** ve standardním WebGL (chybí Storage Buffer Objects) 

---

### **Praktické doporučení**

Pro **fungující implementaci** použijte **Transform Feedback** ve WebGL 2:

```javascript
// Inicializace
const tf = gl.createTransformFeedback();
gl.bindTransformFeedback(gl.TRANSFORM_FEEDBACK, tf);
gl.bindBufferBase(gl.TRANSFORM_FEEDBACK_BUFFER, 0, positionBuffer2);

// Průchod 1: Simulace
gl.useProgram(simulationProgram);
gl.beginTransformFeedback(gl.POINTS);
gl.drawArrays(gl.POINTS, 0, particleCount);
gl.endTransformFeedback();

// Průchod 2: Kreslení
gl.useProgram(renderProgram);
gl.bindBuffer(gl.ARRAY_BUFFER, positionBuffer2);
gl.drawArrays(gl.POINTS, 0, particleCount);
```

---

### **WebGPU jako budoucnost**

Pokud můžete cílit na moderní prohlížeče, **WebGPU** nabízí:
- Nativní compute shadery
- Storage buffery (přímý zápis do bufferů)
- 3D textury s plnou read/write podporou
- **Ale**: Stále experimentalní, ne všechny prohlížeče mají plnou podporu 

**Závěr**: Pro produkční částicové simulace **použijte WebGL 2 + Transform Feedback**. Pro vývoj experimentů s compute shadery zkuste WebGPU.
