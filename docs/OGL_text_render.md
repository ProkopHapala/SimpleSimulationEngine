## User

I want to desing most efficient text rendering system using OpenGL, perhaps using bilboard rendering and UV maping of texture atlas with font. I already have some components (like bilboards, instanced rendering etc. - particle system), however I would like to extend it to render text convenietly.
I also have 1D texture atlas where each letter offset is some multiple of its ASCII code so we can simply use fontSizeX*char_value to calculate texture offset

I was thinking that one option is create bilboard for each letter but that would be perhaps inefficient.
Perhas better is create bilboard for each line (most labels are single-line), so we can map multiple letters on single bilboard (the bilboard will be rectangles with size nlettter*fontSizeX). We then need to compute the coordinates in fragment shader reading from some attribute buffers the char_value for n-th char in the label. Can we do this easily and efficiently?

We cannot use the old fixed OpenGL1 pipeline !!!! we need to reimplement it. What we propably mentioned in the doc that you should look how the setup of texture coordinates from texture atlas should be made by looking how it id done in the lagacy code Draw::drawText function, but what we need to actually do is make some new text rendingin system, where we create quad bilboards. Each letter (char) is simply number, which defines offset in our texture atlas, our texture atlas is 1D. Therefore our text rendering system can easily read directly the letters as an atrubute buffer or someting and then in shader evaluate proper UV coordinates.
Each plabe is simple rectangulat quad with n-letter, and each letter is defined by one index. Can you desing shadr for it?

---

Here is the design and implementation plan.

You are correct to reject the fixed pipeline. Based on your requirements, there are two modern ways to achieve this. Since you asked specifically for the **One Quad Per Label** approach (mapping multiple letters onto a single billboard), I will design the shader for that first.

However, I will also provide the **Instanced Character** approach, because usually, that is actually *more* efficient for GPUs than the "One Quad Per Label" approach due to texture cache locality and shader complexity.

### Approach 1: One Quad Per Label (Your Request)

In this approach, you render **one instance** (one quad) representing the whole word/sentence.
**The Challenge:** You cannot pass a variable-length string as a vertex attribute.
**The Solution:** Use a **TBO (Texture Buffer Object)**. This allows the shader to read the string data (sequence of ASCII codes) from a buffer like a 1D texture.

#### Data Structures
1.  **String Buffer (TBO):** A large buffer containing integers (ASCII codes) for all your labels back-to-back.
2.  **Instance Data:** For each label, you pass:
    *   World Position
    *   `StringOffset` (Where this label starts in the TBO)
    *   `StringLength` (How many characters)

#### The Shader Design

**Vertex Shader:** Pass the geometry and setup the interpolation.
```glsl
#version 330 core

// Instance attributes
layout(location = 0) in vec3 aPos;          // Position of the whole label center/corner
layout(location = 1) in vec2 aSize;         // Size of the label (Width = numChars * charWidth, Height = charHeight)
layout(location = 2) in int  aStringOffset; // Index in TBO where this string starts
layout(location = 3) in int  aStringLength; // Number of chars

// Standard uniforms
uniform mat4 uVP; // ViewProjection
uniform vec3 uCameraRight; // For billboarding
uniform vec3 uCameraUp;    // For billboarding

out vec2 vUV;
flat out int vStringOffset;
flat out int vStringLength;

void main() {
    // Standard Billboard Math (Spherical or Cylindrical)
    // Generating quad vertices locally based on vertexID or input quad geometry
    // Assuming standard 0..1 quad input for simplicity here
    vec3 center = aPos;
    
    // Calculate vertex position relative to center based on UV (0..1)
    // If you use a VBO for the quad, use those attributes. 
    // Here is a procedural quad generation for simplicity:
    vec2 uv = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2); // Generates 0,0 -> 2,0 -> 0,2 covers a quad
    
    vec3 vertexPos = center 
        + uCameraRight * (uv.x - 0.5) * aSize.x 
        + uCameraUp    * (uv.y - 0.5) * aSize.y;

    vUV = uv; // 0..1
    vStringOffset = aStringOffset;
    vStringLength = aStringLength;
    
    gl_Position = uVP * vec4(vertexPos, 1.0);
}
```

**Fragment Shader:** The magic happens here. We use the UV coordinate to find *which* character we are on, fetch that character from the TBO, and calculate the 1D Atlas offset.

```glsl
#version 330 core

in vec2 vUV;
flat in int vStringOffset;
flat in int vStringLength;

// The 1D Font Texture Atlas
uniform sampler2D uFontAtlas;
// The Buffer containing the text strings (ASCII codes)
uniform usamplerBuffer uTextBuffer; 

// Atlas info
uniform float uCharStep; // = 1.0 / TotalCharactersInAtlas (e.g., 1.0 / 128.0)

out vec4 FragColor;

void main() {
    // 1. Determine which character index in the string we are currently rendering
    // vUV.x goes from 0.0 to 1.0 across the whole word.
    float virtualIndex = vUV.x * float(vStringLength);
    
    // The integer part tells us which letter index (0, 1, 2...)
    int charIndex = int(virtualIndex);
    
    // The fractional part tells us the U coordinate INSIDE that specific letter
    float charLocalU = fract(virtualIndex);

    // 2. Fetch the ASCII code from the TBO
    // We assume the buffer holds uints representing ASCII characters
    uint asciiCode = texelFetch(uTextBuffer, vStringOffset + charIndex).r;

    // 3. Map to 1D Texture Atlas
    // Your atlas logic: Offset = ASCII * CharWidth
    // We add charLocalU scaled by the width of one character to map the texture
    float atlasU = (float(asciiCode) * uCharStep) + (charLocalU * uCharStep);
    
    // Sample the font
    // Assuming 1D atlas, V can be 0.5 or ignored depending on texture setup
    vec4 textColor = texture(uFontAtlas, vec2(atlasU, 0.5));

    // Discard transparent pixels (classic text rendering)
    if(textColor.a < 0.1) discard;

    FragColor = textColor; // Add color uniform mixing here if needed
}
```

---

### Approach 2: Instanced Characters (Recommended & Most Efficient)

Since you stated: *"I already have components like instanced rendering - particle system"*.

**This is theoretically the most efficient method.**
Instead of 1 Quad = 1 Label, you do **1 Instance = 1 Character**.
You treat text exactly like a particle system.

*   **Efficiency:** This avoids the dependent texture read (`texelFetch` then `texture`) in the fragment shader seen in Approach 1. The GPU executes the Vertex Shader once per char (very fast) and the Fragment Shader is extremely simple.
*   **Implementation:**
    1.  Generate a buffer of `vec3 position, uint ascii_code` on the CPU for the text label.
    2.  Render the "Particle System" with the Font Atlas texture.

**Vertex Shader (Instanced):**
```glsl
layout(location = 0) in vec3 aCenterPos; // Base position of the label
layout(location = 1) in int  aCharIndex; // Which char in the string (0, 1, 2...)
layout(location = 2) in int  aAsciiCode; // The letter 'A' (65)

uniform vec3 uCameraRight;
uniform vec3 uCameraUp;
uniform float uCharWorldWidth; // Width of char in world space
uniform float uAtlasStep;      // Width of char in texture space (1/128)

out vec2 vUV;

void main() {
    // Calculate world position offset for this specific letter
    vec3 pos = aCenterPos + (uCameraRight * float(aCharIndex) * uCharWorldWidth);
    
    // Billboard logic for the quad corners
    // (Assuming standard Quad vertices 0..1 passed in layout location X)
    vec3 vertexPos = pos 
        + uCameraRight * (QuadPos.x - 0.5) * uCharWorldWidth
        + uCameraUp    * (QuadPos.y - 0.5) * uCharWorldHeight;

    gl_Position = uVP * vec4(vertexPos, 1.0);

    // UV Calculation
    // Base offset based on ASCII + local quad UV scaled down
    float uBase = float(aAsciiCode) * uAtlasStep;
    vUV.x = uBase + (QuadPos.x * uAtlasStep);
    vUV.y = QuadPos.y; // Or 1.0 - QuadPos.y
}
```

---

### Design Document Notes (Correction)

Here is how you should write the note in your design doc:

> **Text Rendering Pipeline Specification:**
> 
> 1.  **Pipeline Type:** Modern Programmable Pipeline (OpenGL 3.3+ / 4.x). Fixed-function (legacy `glBegin/glEnd` or `glTexCoord`) is strictly prohibited.
> 2.  **Texture Source:** A 1D Texture Atlas containing all glyphs.
>     *   Mapping Logic: `UV.u = (ASCII_Code * Glyph_Width_Normalized) + Local_Fragment_U`.
> 3.  **Rendering Technique:**
>     *   **Option A (String-Quads):** Render one billboard quad per string.
>         *   Store string data (ASCII indices) in a **Texture Buffer Object (TBO)**.
>         *   Vertex Shader outputs `StringOffset` and `Length`.
>         *   Fragment Shader calculates the current character index via `floor(UV.x * Length)`, fetches the ASCII code from the TBO using `texelFetch`, and computes the final Atlas UV.
>     *   **Option B (Instanced Glyphs - Preferred):** Utilize the existing Particle System architecture.
>         *   Render one quad instance per character.
>         *   Instance attributes: `WorldPosition`, `ASCII_Code`.
>         *   Vertex Shader computes the texture offset directly.
> 4.  **Geometry:** Quads are strictly View-Aligned (Billboards) using Camera Right/Up vectors in the Vertex Shader.