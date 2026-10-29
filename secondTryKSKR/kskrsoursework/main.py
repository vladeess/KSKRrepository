# main.py
import config
from mesh_generator import generate_mesh
from fem_solver import assemble_global_stiffness, apply_bcs_and_solve
from postprocessing import calculate_stresses_and_export
from postprocessing import visualize_in_python  # Импортируем функцию визуализации


def main():
    print("=== ЭТАП 1: ГЕНЕРАЦИЯ СЕТКИ ===")
    generate_mesh(show_gui=False)

    print("\n=== ЭТАП 2: СБОРКА МАТРИЦЫ ЖЕСТКОСТИ ===")
    K_global, points, tets = assemble_global_stiffness()

    print("\n=== ЭТАП 3: РЕШЕНИЕ СЛАУ ===")
    U_full = apply_bcs_and_solve(K_global, points)

    print("\n=== ЭТАП 4: ПОСТПРОЦЕССИНГ ===")
    calculate_stresses_and_export(points, tets, U_full, output_file="results.vtu")

    print("\n=== ЭТАП 5: ВИЗУАЛИЗАЦИЯ ===")
    # Запускаем 3D окно только после того, как файл results.vtu обновлен!
    visualize_in_python(vtu_file="results.vtu", scale_factor=50.0)


if __name__ == "__main__":
    main()