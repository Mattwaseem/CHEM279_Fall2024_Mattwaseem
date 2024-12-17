// Premptively saves a bunch of variables so that we can use them and reference them throughout
let scene;
let camera;
let renderer;
let controls;

let merged_h_geo;
let O_mat;
let H_mat;

let instance_oxy_m;
let instance_h_m;
let O_color = [1.0, 0.0, 0.0];
let H_color = [1.0, 1.0, 1.0];
let world_scale = 15.0;

let sphereGeometry;
let box_geometry;
let box_material;
let edge_geometry;
let line_geometry;

// Function to save the camera position and orientation to localStorage
function saveCameraState() {
    const cameraState = {
        position: {
            x: camera.position.x,
            y: camera.position.y,
            z: camera.position.z
        },
        rotation: {
            x: camera.rotation.x,
            y: camera.rotation.y,
            z: camera.rotation.z
        }
    };
    localStorage.setItem('cameraState', JSON.stringify(cameraState));  // Save to localStorage
}

// Function to restore the camera position and orientation from localStorage
function restoreCameraState() {
    const cameraState = localStorage.getItem('cameraState');
    if (cameraState) {
        const { position, rotation } = JSON.parse(cameraState);
        camera.position.set(position.x, position.y, position.z);
        camera.rotation.set(rotation.x, rotation.y, rotation.z);
    }
}

// Initializes threeJS
function initThreeJS() {

    // Setup scene and camera
    scene = new THREE.Scene();
    init_camera();

    init_box();

    init_spheres();

    // OrbitControls
    controls = new THREE.OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.25;
    controls.screenSpacePanning = false;

    addLighting();

    // Start the animation loop
    animate();
}

function init_camera() {
    camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.001, 10000);
    restoreCameraState();
    renderer = new THREE.WebGLRenderer();

    // Optionally, set a default camera position if no state was saved
    if (!localStorage.getItem('cameraState')) {
        camera.position.z = 5;  // Default position if no previous state is saved
    }

    // Fill window of threejs container
    renderer.setSize(window.innerWidth, window.innerHeight);
    document.getElementById('threejs-container').appendChild(renderer.domElement);
}

function init_spheres() {

    O_mat = new THREE.RawShaderMaterial({
        uniforms: {
            'map': { value: new THREE.TextureLoader().load('textures/sprites/circle.png') },
            'part_color': { value: O_color }
        },
        vertexShader: document.getElementById('vshader').textContent,
        fragmentShader: document.getElementById('fshader').textContent,
        depthTest: true,
        depthWrite: true
    });

    H_mat = new THREE.RawShaderMaterial({
        uniforms: {
            'map': { value: new THREE.TextureLoader().load('textures/sprites/circle.png') },
            'part_color': { value: H_color }
        },
        vertexShader: document.getElementById('vshader').textContent,
        fragmentShader: document.getElementById('fshader').textContent,
        depthTest: true,
        depthWrite: true
    });

    let oxy_geo = new THREE.SphereGeometry(0.05, 8, 8);
    let h1_geo = new THREE.SphereGeometry(0.03, 8, 8);
    let h2_geo = new THREE.SphereGeometry(0.03, 8, 8);

    let rads = ((109.5/2.0 - 90.0)*Math.PI/180.0);
    h1_geo.translate(-1.0 / world_scale * Math.cos(rads), 1.0 / world_scale * Math.sin(rads), 0.0);
    h2_geo.translate(1.0 / world_scale * Math.cos(rads), 1.0 / world_scale * Math.sin(rads), 0.0);
    merged_h_geo = THREE.BufferGeometryUtils.mergeBufferGeometries([
        h1_geo, h2_geo
    ]);

    instance_oxy_m = new THREE.InstancedMesh(oxy_geo, O_mat, 0);
    instance_h_m = new THREE.InstancedMesh(merged_h_geo, H_mat, 0);

    scene.add(instance_oxy_m);
    scene.add(instance_h_m);
}

function init_box() {

    let w = 2;
    let h = 2;
    let d = 2;
    // // Create the box geometry
    box_geometry = new THREE.BoxGeometry(w, h, d);  // Size of the box
    box_material = new THREE.LineBasicMaterial({ color: 0x00ff00 });  // Color for the edges

    // Use EdgesGeometry to get the edges
    edge_geometry = new THREE.EdgesGeometry(box_geometry);

    // Line segments to be rendered
    line_geometry = new THREE.LineSegments(edge_geometry, box_material);
    const quaternion = new THREE.Quaternion(0.0, 0.0, 0.0, 1.0);
    // const quaternion = new THREE.Quaternion(-0.92388, -0.22068, 0.22068, 0.22068);
    line_geometry.rotation.setFromQuaternion(quaternion);
    scene.add(line_geometry);
    // update_box_size(w, h, d);
}

function addLighting() {
    // 1. Ambient Light (soft, uniform light)
    const ambientLight = new THREE.AmbientLight(0x0000ff, 1);  // Color and intensity
    scene.add(ambientLight);

    // 2. Directional Light (hard, directional light)
    const directionalLight = new THREE.DirectionalLight(0xccaa00, 1);  // Color and intensity
    directionalLight.position.set(0, 10, 0);  // Set the light position
    scene.add(directionalLight);

    // Pass lighting information as uniforms to the material
    O_mat.uniforms.dirLightColor = { value: directionalLight.color };
    O_mat.uniforms.dirLightDir = { value: directionalLight.position.clone().normalize() };  // Normalize the direction
    O_mat.uniforms.dirLightPos = { value: directionalLight.position.clone() };  // Normalize the direction
    O_mat.uniforms.ambientLightColor = { value: ambientLight.color };
    H_mat.uniforms.dirLightColor = { value: directionalLight.color };
    H_mat.uniforms.dirLightDir = { value: directionalLight.position.clone().normalize() };  // Normalize the direction
    H_mat.uniforms.dirLightPos = { value: directionalLight.position.clone() };  // Normalize the direction
    H_mat.uniforms.ambientLightColor = { value: ambientLight.color };
}

let sim_running = false;
function animate() {
    // Used by nodejs, and WebGL. It's so get the rendering engine ready to render the next frame.
    requestAnimationFrame(animate);
    const time = performance.now() * 0.0005;
    // O_color,at.uniforms['time'].value = time;

    // // Step the sim once we've started it.
    if (sim_running)
        step_sim();

    // Update OrbitControls & save current camera position to local storage so we can
    // refresh the page as many times as we like without needing to reset the camera.
    controls.update();
    saveCameraState();

    // Render the scene
    renderer.render(scene, camera);
}

// Handle window resizing
window.addEventListener('resize', () => {
    renderer.setSize(window.innerWidth, window.innerHeight);
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
});

// https://threejs.org/docs/scenes/geometry-browser.html#SphereGeometry
function update_point_size(size) {
    if (instance_oxy_m) {
        scene.remove(instance_oxy_m);
        instance_oxy_m.geometry.dispose();
    }

    if (instance_h_m) {
        scene.remove(instance_h_m);
        instance_h_m.geometry.dispose();
    }
    let oxy_geo = new THREE.SphereGeometry(size, 8, 8);
    let h1_geo = new THREE.SphereGeometry(size * 0.6, 8, 8);
    let h2_geo = new THREE.SphereGeometry(size * 0.6, 8, 8);
    h1_geo.translate(size * 0.9, -size * 0.5, 0.0);
    h2_geo.translate(-size * 0.9, -size * 0.5, 0.0);
    merged_h_geo = THREE.BufferGeometryUtils.mergeBufferGeometries([
        h1_geo, h2_geo
    ]);

    instance_oxy_m = new THREE.InstancedMesh(oxy_geo, O_mat, parameter_container.num_points);
    instance_h_m = new THREE.InstancedMesh(merged_h_geo, H_mat, parameter_container.num_points);

    // Set the position attribute for instancing
    instance_oxy_m.geometry.setAttribute('instancePosition', new THREE.InstancedBufferAttribute(curr_positions, 3));
    instance_h_m.geometry.setAttribute('instancePosition', new THREE.InstancedBufferAttribute(curr_positions, 3));

    scene.add(instance_oxy_m);
    scene.add(instance_h_m);
}

function update_box_size(num_parts) {

    let part_diff = 4.0;
    let root = Math.cbrt(num_parts);
    let side_length = root;
    if (root * root * root != num_parts)
        side_length += 1;
    let world_size = (side_length + 1) * part_diff * inv_world_scale;
    let center = world_size * 0.5 - part_diff;
    
    box_geometry.dispose();
    box_geometry = new THREE.BoxGeometry(world_size, world_size, world_size);
    edge_geometry = new THREE.EdgesGeometry(box_geometry);
    line_geometry.geometry = edge_geometry;
    const quaternion = new THREE.Quaternion(0.0, 0.0, 0.0, 1.0);
    // const quaternion = new THREE.Quaternion(-0.92388, -0.22068, 0.22068, 0.22068);
    line_geometry.rotation.setFromQuaternion(quaternion);

    // Used to live in Sim file, but this will have to do.
    fetch(`/set/sim/box?w=${world_size}&h=${world_size}&d=${world_size}`)
        .then(response => response.json())// Print errors if there are any.
        .catch(error => {
            console.error('Error:', error);
            document.getElementById("error").textContent = "Error updating box parameters";
        });
}

function update_part_cols(o_col, h_col) {

    O_color = [o_col[0] / 255.0, o_col[1] / 255.0, o_col[2] / 255.0];
    H_color = [h_col[0] / 255.0, h_col[1] / 255.0, h_col[2] / 255.0];

    O_mat.uniforms['part_color'].value = O_color;
    H_mat.uniforms['part_color'].value = H_color;
}

function update_lighting_cols() {
    let dir_col = parameter_container.dir_light_col;
    let amb_col = parameter_container.ambient_light_col;

    O_mat.uniforms['dirLightDir'].value = parameter_container.dir_light_dir.normalize();
    O_mat.uniforms['dirLightPos'].value = parameter_container.dir_light_pos;
    O_mat.uniforms['dirLightColor'].value = [dir_col[0] / 255.0, dir_col[1] / 255.0, dir_col[2] / 255.0];
    O_mat.uniforms['ambientLightColor'].value = [amb_col[0] / 255.0, amb_col[1] / 255.0, amb_col[2] / 255.0];

    H_mat.uniforms['dirLightDir'].value = parameter_container.dir_light_dir.normalize().toArray();
    H_mat.uniforms['dirLightPos'].value = parameter_container.dir_light_pos.toArray();
    H_mat.uniforms['dirLightColor'].value = [dir_col[0] / 255.0, dir_col[1] / 255.0, dir_col[2] / 255.0];
    H_mat.uniforms['ambientLightColor'].value = [amb_col[0] / 255.0, amb_col[1] / 255.0, amb_col[2] / 255.0];
}


initThreeJS(); // Always keep at bottom