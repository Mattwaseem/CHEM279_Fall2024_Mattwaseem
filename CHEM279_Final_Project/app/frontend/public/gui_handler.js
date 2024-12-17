
// https://github.com/dataarts/dat.gui/blob/master/example.html
let parameter_container = {
    label: 'Water Simulation',
    sim_scale: 0.001,
    point_size: 0.02,
    height: 10,
    width: 10,
    depth: 10,
    num_points: 300,
    init_temp: 100,
    H_color: [0, 0, 0], // RGB array
    O_color: [255, 0, 0], // RGB array
    ambient_light_col: [0, 0, 0], // CSS string
    dir_light_col: [0, 128, 255], // RGB
    dir_light_pos: new THREE.Vector3(10, 0, 0), // XYZ
    dir_light_dir: new THREE.Vector3(1, 0, 0), // XYZ (gets normalized)
};

var gui = new dat.gui.GUI();
gui.remember(parameter_container);
// elem = gui.add(sim_name, 'label');
// elem.domElement.querySelector('input').disabled = true;

gui.add(parameter_container, 'point_size')
    .min(0.01)
    .max(0.1)
    .onChange(function (new_value) {
        update_point_size(new_value);
});
var spawn_folder = gui.addFolder("Spawner");

// spawn_folder.add(parameter_container, "sim_scale").min(0.000001).max(0.01);

// let box_dims_folder = spawn_folder.addFolder("Box Dimensions");
// box_dims_folder.add(parameter_container, 'width').min(0.01).max(10)
//     .onChange(function (new_value) {
//         // update_box_size(new_value, parameter_container.height, parameter_container.depth);
//     });
// box_dims_folder.add(parameter_container, 'height').min(0.01).max(10)
//     .onChange(function (new_value) {
//         // update_box_size(parameter_container.width, new_value, parameter_container.depth);
//     });
// box_dims_folder.add(parameter_container, 'depth').min(0.01).max(10)
//     .onChange(function (new_value) {
//         // update_box_size(parameter_container.width, parameter_container.height, new_value);
//     });

spawn_folder.add(parameter_container, 'num_points').min(1).max(500).step(1)
    .onChange(function (new_val) {
        initSimulation(parameter_container.num_points, parameter_container.init_temp);
        update_box_size(new_val);
    });
spawn_folder.add(parameter_container, 'init_temp').min(0).max(500)
    .onChange(function (new_val) {
        initSimulation(parameter_container.num_points, parameter_container.init_temp);
    });

var color_folder = gui.addFolder('Colors');
color_folder.addColor(parameter_container, 'O_color')
    .onChange(function (new_value) {
        update_part_cols(new_value, parameter_container.H_color);
    });
color_folder.addColor(parameter_container, 'H_color')
    .onChange(function (new_value) {
        update_part_cols(parameter_container.O_color, new_value);
    });

var light_folder = gui.addFolder('Lighting');
light_folder.addColor(parameter_container, 'ambient_light_col')
    .onChange(update_lighting_cols);
light_folder.addColor(parameter_container, 'dir_light_col')
    .onChange(update_lighting_cols);

var dir_pos_folder = light_folder.addFolder("Directional Light Position");
dir_pos_folder.add(parameter_container.dir_light_pos, 'x', -10, 10)
    .onChange(update_lighting_cols);
dir_pos_folder.add(parameter_container.dir_light_pos, 'y', -10, 10)
    .onChange(update_lighting_cols);
dir_pos_folder.add(parameter_container.dir_light_pos, 'z', -10, 10)
    .onChange(update_lighting_cols);

var dir_dir_folder = light_folder.addFolder("Directional Light Direction");
dir_dir_folder.add(parameter_container.dir_light_dir, 'y', -1, 1)
    .onChange(update_lighting_cols);
dir_dir_folder.add(parameter_container.dir_light_dir, 'x', -1, 1)
    .onChange(update_lighting_cols);
dir_dir_folder.add(parameter_container.dir_light_dir, 'z', -1, 1)
    .onChange(update_lighting_cols);

parameter_container['Step Simulation'] = function () {
    step_sim();
};
parameter_container['Start Simulation'] = function () {
    start_simulation();
};
parameter_container['Stop Simulation'] = function () {
    stop_simulation();
};

gui.add(parameter_container, 'Step Simulation');
gui.add(parameter_container, 'Start Simulation');
gui.add(parameter_container, 'Stop Simulation');

function disable_spawner_gui() {
    
}