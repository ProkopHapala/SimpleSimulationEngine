
export class CatalogItem {
    constructor(id, name, kind) {
        this.id = id;
        this.name = name;
        this.kind = kind;
    }
}

export class Material extends CatalogItem {
    constructor(id, name, density, Spull, Spush, Kpull, Kpush, reflectivity, Tmelt) {
        super(id, name, 0); 
        this.density = density;      // [kg/m3]
        this.Spull = Spull;          // [Pa] Strength
        this.Spush = Spush;          // [Pa]
        this.Kpull = Kpull;          // [Pa] elastic modulus
        this.Kpush = Kpush;          // [Pa]
        this.reflectivity = reflectivity;
        this.Tmelt = Tmelt;          // [K]
    }
}

export class StickMaterial extends CatalogItem {
    constructor(id, name, materialId, diameter, wallThickness) {
        super(id, name, 1);
        this.materialId = materialId;
        this.diameter = diameter;
        this.wallThickness = wallThickness;
        
        // Physical parameters (updated from Material)
        this.area = 0;
        this.linearDensity = 0;
        this.Kpull = 0;
        this.Kpush = 0;
        this.Spull = 0;
        this.Spush = 0;
        this.reflectivity = 0;
        this.Tmelt = 0;
        this.damping = 0.1; // Default damping
        this.preStrain = 0.0;
    }

    update(materials) {
        const mat = materials[this.materialId];
        if (!mat) return;
        const din = this.diameter - 2 * this.wallThickness;
        this.area = Math.PI * (this.diameter * this.diameter - din * din) * 0.25;
        this.linearDensity = this.area * mat.density;
        this.Kpull = this.area * mat.Kpull;
        this.Kpush = this.area * mat.Kpush;
        this.Spull = this.area * mat.Spull;
        this.Spush = this.area * mat.Spush;
        this.reflectivity = mat.reflectivity;
        this.Tmelt = mat.Tmelt;
    }
}

export class PanelLayer {
    constructor(materialId, thickness) {
        this.materialId = materialId;
        this.thickness = thickness;
    }
}

export class PanelMaterial extends CatalogItem {
    constructor(id, name, stickMaterialId) {
        super(id, name, 2);
        this.layers = [];
        this.areaDensity = 0;
        this.stickMaterialId = stickMaterialId;
    }

    addLayer(materialId, thickness) {
        this.layers.push(new PanelLayer(materialId, thickness));
    }

    evalAreaDensity(materials) {
        this.areaDensity = 0;
        for (const layer of this.layers) {
            const mat = materials[layer.materialId];
            if (mat) {
                this.areaDensity += layer.thickness * mat.density;
            }
        }
    }
}

export class SpaceCraftWorkshop {
    constructor() {
        this.materials      = [];
        this.stickMaterials = [];
        this.panelMaterials = [];
        this.materialsMap   = new Map();
        this.stickMaterialsMap = new Map();
        this.initDefaults();
    }

    initDefaults() {
        this.addMaterial("Steel",       7800.0, 400e6,  400e6,  200e9, 200e9, 0.5, 1800.0 );
        this.addMaterial("Aluminum",    2700.0, 200e6,  200e6,  70e9,  70e9,  0.8, 900.0  );
        this.addMaterial("CarbonFiber", 1600.0, 2000e6, 2000e6, 150e9, 150e9, 0.1, 3000.0 );
        this.addStickMaterial("SteelGirder",       "Steel",       0.1,  0.005 );
        this.addStickMaterial("AluminumGirder",    "Aluminum",    0.1,  0.005 );
        this.addStickMaterial("CarbonFiberGirder", "CarbonFiber", 0.1,  0.005 );
        this.addStickMaterial("rope",              "CarbonFiber", 0.05, 0.05  ); // Solid rope
    }

    addMaterial(name, density, Spull, Spush, Kpull, Kpush, reflectivity, Tmelt) {
        const id = this.materials.length;
        const mat = new Material(id, name, density, Spull, Spush, Kpull, Kpush, reflectivity, Tmelt);
        this.materials.push(mat);
        this.materialsMap.set(name, id);
        return id;
    }

    addStickMaterial(name, matName, diameter, wallThickness) {
        const matId = this.materialsMap.get(matName);
        if (matId === undefined) throw new Error(`Material ${matName} not found`);
        const id = this.stickMaterials.length;
        const stick = new StickMaterial(id, name, matId, diameter, wallThickness);
        stick.update(this.materials);
        this.stickMaterials.push(stick);
        this.stickMaterialsMap.set(name, id);
        return id;
    }

    addPanelMaterial(name, stickMatName, layers = []) {
        const stickId = this.stickMaterialsMap.get(stickMatName);
        const id = this.panelMaterials.length;
        const panel = new PanelMaterial(id, name, stickId);
        for (const l of layers) {
            const matId = this.materialsMap.get(l.name);
            panel.addLayer(matId, l.thickness);
        }
        panel.evalAreaDensity(this.materials);
        this.panelMaterials.push(panel);
        return id;
    }
}
