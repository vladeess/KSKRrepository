import numpy as np
import meshio
import pyvista as pv
import config
from fem_solver import get_elasticity_matrix, calculate_tet_volume_and_B


def calculate_stresses_and_export(points, tets, U_full, output_file="results.vtu"):
    print("[1/3] Вычисление напряжений и деформаций для каждого элемента...")
    num_elements = len(tets)

    stresses = np.zeros((num_elements, 6))
    von_mises = np.zeros(num_elements)

    D_outer = get_elasticity_matrix(config.MATERIALS["outer_skin"]["E"], config.MATERIALS["outer_skin"]["nu"])
    D_core = get_elasticity_matrix(config.MATERIALS["core"]["E"], config.MATERIALS["core"]["nu"])
    D_inner = get_elasticity_matrix(config.MATERIALS["inner_skin"]["E"], config.MATERIALS["inner_skin"]["nu"])

    for i, tet in enumerate(tets):
        tet_nodes = points[tet]
        z_center = np.mean(tet_nodes[:, 2])
        if z_center <= config.T_OUTER:
            D = D_outer
        elif z_center <= (config.T_OUTER + config.T_CORE):
            D = D_core
        else:
            D = D_inner

        _, B = calculate_tet_volume_and_B(tet_nodes)

        u_e = np.zeros(12)
        for j in range(4):
            u_e[j * 3] = U_full[tet[j] * 3]
            u_e[j * 3 + 1] = U_full[tet[j] * 3 + 1]
            u_e[j * 3 + 2] = U_full[tet[j] * 3 + 2]

        eps = B @ u_e
        sig = D @ eps
        stresses[i, :] = sig

        sx, sy, sz, sxy, syz, szx = sig
        vm = np.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2 + 6 * (sxy ** 2 + syz ** 2 + szx ** 2)))
        von_mises[i] = vm

    print("\n[2/3] Формирование данных для отчета...")
    max_vm = np.max(von_mises)
    max_vm_element = np.argmax(von_mises)

    sig_max = stresses[max_vm_element]
    nodes_max = tets[max_vm_element]

    print("\n================ РЕЗУЛЬТАТЫ ДЛЯ ОТЧЕТА ================")
    print(f"Самый нагруженный элемент: №{max_vm_element} (состоит из узлов {nodes_max})")
    print(f"Эквивалентное напряжение Мизеса: {max_vm / 1e6:.3f} МПа")
    print("\nТензор напряжений {Sigma} (Па):")
    print(f"  Sig_x  = {sig_max[0]:11.3E}")
    print(f"  Sig_y  = {sig_max[1]:11.3E}")
    print(f"  Sig_z  = {sig_max[2]:11.3E}")
    print(f"  Tau_xy = {sig_max[3]:11.3E}")
    print(f"  Tau_yz = {sig_max[4]:11.3E}")
    print(f"  Tau_zx = {sig_max[5]:11.3E}")

    displacements_3d = U_full.reshape((-1, 3))
    disp_magnitudes = np.linalg.norm(displacements_3d, axis=1)
    max_node_idx = np.argmax(disp_magnitudes)

    print("\nВекторы узловых перемещений (выборка):")
    print("  Узел |    Ux (мм) |    Uy (мм) |    Uz (мм) | Описание")
    print("  --------------------------------------------------------")

    fixed_printed = 0
    for i in range(len(displacements_3d)):
        if np.allclose(displacements_3d[i], 0.0, atol=1e-10):
            ux, uy, uz = displacements_3d[i] * 1000
            print(f"  {i:4d} | {ux:10.5f} | {uy:10.5f} | {uz:10.5f} | Заделка (Y=0)")
            fixed_printed += 1
            if fixed_printed == 2:
                break

    ux, uy, uz = displacements_3d[max_node_idx] * 1000
    print(f"  {max_node_idx:4d} | {ux:10.5f} | {uy:10.5f} | {uz:10.5f} | МАКСИМАЛЬНЫЙ ПРОГИБ")
    print("=======================================================\n")

    cells = [("tetra", tets)]
    mesh = meshio.Mesh(
        points,
        cells,
        point_data={"Displacement": displacements_3d},
        cell_data={
            "Sigma_x": [stresses[:, 0]],
            "Sigma_y": [stresses[:, 1]],
            "Sigma_z": [stresses[:, 2]],
            "Von_Mises": [von_mises]
        }
    )

    print(f"[3/3] Сохранение результатов в файл {output_file}...")
    mesh.write(output_file)
    print("Готово!")


def visualize_in_python(vtu_file="results.vtu", scale_factor=50.0):
    mesh = pv.read(vtu_file)
    warped_mesh = mesh.warp_by_vector("Displacement", factor=scale_factor)

    plotter = pv.Plotter(title="МКЭ: Напряженно-деформированное состояние")
    plotter.add_mesh(
        warped_mesh,
        scalars="Von_Mises",
        cmap="jet",
        show_edges=True,
        edge_color="white",
        line_width=0.5
    )

    plotter.add_axes()
    print("Запуск 3D-окна... Используйте мышь для вращения детали.")
    plotter.show()