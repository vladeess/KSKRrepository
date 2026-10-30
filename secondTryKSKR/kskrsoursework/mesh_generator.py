import gmsh
import math


def generate_mesh(output_file="mesh.msh"):
    """
    Генератор 3D-сетки для многослойной сэндвич-панели с вырезами и усилителями.
    Используется библиотека Gmsh (OpenCASCADE ядро).
    """
    gmsh.initialize()
    gmsh.model.add("sandwich_panel")

    # ==========================================
    # 1. ГЕОМЕТРИЧЕСКИЕ ПАРАМЕТРЫ ПАНЕЛИ
    # ==========================================
    L = 6.0  # Длина панели (по оси X)
    H = 3.0  # Высота панели (по оси Y)

    t_out = 0.05  # Толщина наружной стальной обшивки (по Z)
    t_core = 0.20  # Толщина внутреннего полимерного утеплителя
    t_in = 0.05  # Толщина внутренней стальной обшивки

    # ==========================================
    # 2. ФОРМИРОВАНИЕ БАЗОВЫХ СЛОЕВ (Объемы)
    # ==========================================
    # gmsh.model.occ.addBox(x, y, z, dx, dy, dz)
    v_out = gmsh.model.occ.addBox(0, 0, 0, L, H, t_out)
    v_core = gmsh.model.occ.addBox(0, 0, t_out, L, H, t_core)
    v_in = gmsh.model.occ.addBox(0, 0, t_out + t_core, L, H, t_in)

    # ==========================================
    # 3. ДОБАВЛЕНИЕ УСИЛИТЕЛЕЙ (Ребра жесткости)
    # ==========================================
    z_stiff = t_out + t_core + t_in  # Координата начала усилителей
    t_stiff = 0.10  # Вылет усилителя (толщина)
    w_stiff = 0.15  # Ширина ребра жесткости

    # Вертикальные ребра жесткости на внутренней обшивке
    s1 = gmsh.model.occ.addBox(0, 0, z_stiff, w_stiff, H, t_stiff)
    s2 = gmsh.model.occ.addBox(2.0, 0, z_stiff, w_stiff, H, t_stiff)
    s3 = gmsh.model.occ.addBox(4.0, 0, z_stiff, w_stiff, H, t_stiff)
    s4 = gmsh.model.occ.addBox(L - w_stiff, 0, z_stiff, w_stiff, H, t_stiff)

    # ==========================================
    # 4. ФОРМИРОВАНИЕ ВЫРЕЗОВ (Булевы инструменты)
    # ==========================================
    # Задаем глубину инструментов с запасом, чтобы они гарантированно прорезали все слои
    cut_z = -0.5
    cut_dz = 2.0

    # Левое окно
    hole1 = gmsh.model.occ.addBox(0.4, 1.0, cut_z, 1.2, 1.2, cut_dz)

    # Центральный вырез (ИСПРАВЛЕНО: теперь это окно с цельным основанием)
    # Основание поднято на y = 0.5 (нижняя кромка не разрезается)
    hole2 = gmsh.model.occ.addBox(2.4, 0.5, cut_z, 1.2, 1.7, cut_dz)

    # Правое окно
    hole3 = gmsh.model.occ.addBox(4.4, 1.0, cut_z, 1.2, 1.2, cut_dz)

    # ==========================================
    # 5. БУЛЕВЫ ОПЕРАЦИИ (Вычитание и Сшивка)
    # ==========================================
    # Собираем все детали конструкции в один список
    parts = [(3, v_out), (3, v_core), (3, v_in), (3, s1), (3, s2), (3, s3), (3, s4)]
    holes = [(3, hole1), (3, hole2), (3, hole3)]

    # Шаг 5.1: Вычитаем инструменты (отверстия) из базовых деталей
    cut_res, _ = gmsh.model.occ.cut(parts, holes, removeObject=True, removeTool=True)

    # Шаг 5.2: Фрагментация (Топологическая сшивка)
    # КРИТИЧЕСКИ ВАЖНО: fragment гарантирует слияние узлов на границах раздела
    # обшивки и утеплителя, исключая расслоение модели при изгибе в МКЭ.
    gmsh.model.occ.fragment(cut_res, [], removeObject=True, removeTool=True)

    gmsh.model.occ.synchronize()

    # ==========================================
    # 6. НАСТРОЙКИ ГЕНЕРАТОРА СЕТКИ
    # ==========================================
    # Параметры плотности сетки (размеры тетраэдров)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.05)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.25)

    # Включаем алгоритмы оптимизации формы элементов (улучшает определитель Якоби)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)

    print("[1/2] Запущена генерация 3D тетраэдральной сетки...")
    gmsh.model.mesh.generate(3)

    print(f"[2/2] Сохранение топологии в файл {output_file}...")
    gmsh.write(output_file)

    gmsh.finalize()
    print("Генерация сетки успешно завершена! Массив узлов сшит корректно.")


if __name__ == "__main__":
    generate_mesh("mesh.msh")