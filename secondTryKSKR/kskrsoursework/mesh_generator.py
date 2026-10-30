import gmsh
import math
import os
import sys

# ИМЯ ФАЙЛА СЕТКИ (Должно совпадать с тем, что ищет ваш интерфейс)
MESH_FILENAME = "panel_mesh.msh"


def generate_mesh(output_file=MESH_FILENAME):
    # 1. ЖЕСТКАЯ ПРОВЕРКА НА БЛОКИРОВКУ ФАЙЛА
    if os.path.exists(output_file):
        try:
            os.remove(output_file)
            print(f"Старый файл '{output_file}' успешно удален.")
        except PermissionError:
            print(f"\n[КРИТИЧЕСКАЯ ОШИБКА] Файл '{output_file}' заблокирован!")
            print("У вас открыто окно интерфейса (main_gui.py) или файл используется!")
            print("ЗАКРОЙТЕ ВСЕ ОКНА с 3D-моделью и запустите генератор снова!")
            sys.exit(1)  # Жестко останавливаем скрипт, чтобы не вводить в заблуждение

    gmsh.initialize()
    gmsh.model.add("sandwich_panel")

    # ==========================================
    # 1. ГЕОМЕТРИЯ ПАНЕЛИ
    # ==========================================
    L = 6.0
    H = 3.0
    t_out = 0.05
    t_core = 0.20
    t_in = 0.05

    v_out = gmsh.model.occ.addBox(0, 0, 0, L, H, t_out)
    v_core = gmsh.model.occ.addBox(0, 0, t_out, L, H, t_core)
    v_in = gmsh.model.occ.addBox(0, 0, t_out + t_core, L, H, t_in)

    z_stiff = t_out + t_core + t_in
    t_stiff = 0.10
    w_stiff = 0.15

    s1 = gmsh.model.occ.addBox(0, 0, z_stiff, w_stiff, H, t_stiff)
    s2 = gmsh.model.occ.addBox(2.0, 0, z_stiff, w_stiff, H, t_stiff)
    s3 = gmsh.model.occ.addBox(4.0, 0, z_stiff, w_stiff, H, t_stiff)
    s4 = gmsh.model.occ.addBox(L - w_stiff, 0, z_stiff, w_stiff, H, t_stiff)

    # ==========================================
    # 2. ФОРМИРОВАНИЕ ВЫРЕЗОВ С ПОРОГОМ
    # ==========================================
    cut_z = -1.0
    cut_dz = 3.0

    hole1 = gmsh.model.occ.addBox(0.4, 1.0, cut_z, 1.2, 1.2, cut_dz)

    # ЦЕНТРАЛЬНОЕ ОКНО С ПОРОГОМ:
    # y = 0.4 означает, что вырез начнется в 40 сантиметрах от пола (сплошной порог)
    # 1.8 - высота самого окна (до уровня 2.2 метра, как и у боковых окон)
    hole2 = gmsh.model.occ.addBox(2.4, 0.4, cut_z, 1.2, 1.8, cut_dz)

    hole3 = gmsh.model.occ.addBox(4.4, 1.0, cut_z, 1.2, 1.2, cut_dz)

    parts = [(3, v_out), (3, v_core), (3, v_in), (3, s1), (3, s2), (3, s3), (3, s4)]
    holes = [(3, hole1), (3, hole2), (3, hole3)]

    cut_res, _ = gmsh.model.occ.cut(parts, holes, removeObject=True, removeTool=True)
    gmsh.model.occ.fragment(cut_res, [], removeObject=True, removeTool=True)
    gmsh.model.occ.synchronize()

    # ==========================================
    # 3. ГЕНЕРАЦИЯ СЕТКИ
    # ==========================================
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.05)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.25)

    print("Генерация 3D сетки (подождите)...")
    gmsh.model.mesh.generate(3)
    gmsh.write(output_file)
    gmsh.finalize()

    print(f"\nГОТОВО! Новая геометрия с порогом сохранена в файл '{output_file}'.")
    print("Теперь можете запускать main_gui.py!")


if __name__ == "__main__":
    generate_mesh()