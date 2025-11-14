// ===================== 
//   THREE JS MAIN
// =====================

// ----  Globals

var renderer;
var scene;
var camera;

var uniforms = {};
var control;
var mousePosOld;

var basicShader;
var mesh;

var str_PixelMap = "vec4( (c_diffuse + c_specular*mat.gloss)*mat.color + vec3(0.1,0.1,0.2)*mat.color, 1.0 );";

var primitivesSrc = "";
var rayTracerTemplate = "";
var shaderFilesReady = Promise.all([
    fetch("Primitives.glslf").then(function(response){ return response.text(); }).then(function(text){ primitivesSrc = text; }),
    fetch("RayTracer.glslf").then(function(response){ return response.text(); }).then(function(text){ rayTracerTemplate = text; })
]);

// ---- Functions

var camQuat = new  THREE.Quaternion();    
function handleMouseMove(event) {
    var dot, eventDoc, doc, body, pageX, pageY;        
    if (mousePosOld) {
        var dx = (event.clientX-mousePosOld.x)*1.0;
        var dy = (event.clientY-mousePosOld.y)*1.0;
        var q = new THREE.Quaternion( -dy*0.002, 0.0, dx*0.002, 1.0 );
        camQuat.multiply(q).normalize();            
    }else{
        mousePosOld = new THREE.Vector3();
    }
    mousePosOld.x = event.clientX;
    mousePosOld.y = event.clientY;
}

function render() {        
    var camMat4 = new THREE.Matrix4(); camMat4.compose ( new THREE.Vector3(0.0,0.0,0.0), camQuat, new THREE.Vector3(1.0,1.0,1.0) );
    var camMat_ = new THREE.Matrix3(); camMat_.getNormalMatrix ( camMat4 );
    //console.log( camMat_.elements +" "+ camQuat.x+" "+ camQuat.y+" "+ camQuat.z+" "+ camQuat.w  );
    uniforms.camMat.value = camMat_;
    
    renderer.render(scene, camera);
    uniforms.time.value += 0.05;
    requestAnimationFrame(render);
}

function selectPrimitive(element){
    console.log(element.value);
    txtScene = "vec4 scene( Ray ray ){\n\tvec4 hit = vec4( POSITIVE_INF, vec3(0.0) );\n\t"+element.value+"\n\tvec2 ts1;\n\tADD( SURF1 );\n\treturn hit;\n}";
    document.getElementById("txtScene").value = txtScene;
    updateShader(document.getElementById("txtScene"));
}

function selectPixelMap(element){
    console.log(element.value);
    str_PixelMap = element.value;
    updateShader(document.getElementById("txtScene"));
}

function updateShader(element){
	//console.log("updateShader");
	//console.log(element.value);

	var shader_code="";
	shader_code += primitivesSrc;
	shader_code += element.value;
	shader_code += rayTracerTemplate.replace("OUTPUT_PIXEL",str_PixelMap);
			
    var material = new THREE.ShaderMaterial({
        uniforms: uniforms,
        vertexShader: basicShader.vertexShader,
        fragmentShader: shader_code,
    });
    		
	try {
	    var mesh = scene.getObjectByName( "Mesh1" );
        scene.remove( mesh );
    }catch(err) {
        console.log( " updateShader cannot remove Mesh1" );
    }; 

	mesh = new THREE.Mesh( new THREE.PlaneBufferGeometry( 100, 100 ), material );
	mesh.name = "Mesh1";
	scene.add( mesh );
}

function init_GLSLScreen(screenBoxId,shaderBoxId){
    screenBox = document.getElementById(screenBoxId);
    document.onmousemove = handleMouseMove;
    console.log(screenBoxId,shaderBoxId );
    scene    = new THREE.Scene();
    camera   = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.1, 1000);
    renderer = new THREE.WebGLRenderer();
    renderer.setClearColor(0x000000, 1.0);
    //renderer.setSize(window.innerWidth, window.innerHeight);
    renderer.setSize( screenBox.clientWidth, screenBox.clientHeight );

    // SHADER
    basicShader = THREE.ShaderLib['normal'];
	uniforms = {
			time      : { value: 1.0 },
			resolution: { value: new THREE.Vector2() },
			camMat    : { value: new THREE.Matrix3() }
	};
	uniforms.resolution.value.x = renderer.domElement.width;
	uniforms.resolution.value.y = renderer.domElement.height;
	shaderFilesReady.then(function(){
	    updateShader( document.getElementById(shaderBoxId) );
	});
			
    camera.position.x = 0.0;
    camera.position.y = 0.0;
    camera.position.z = 100.0;
    camera.lookAt(scene.position);
    
    //document.body.appendChild(renderer.domElement);

    screenBox.appendChild( renderer.domElement );
    control = new function () {
        this.rotationSpeed = 0.005;
        this.scale = 1;
    };
    //addControls(control);
    // call the render function
    render();
}

//window.onload = init; // calls the init function when the window is done loading.
