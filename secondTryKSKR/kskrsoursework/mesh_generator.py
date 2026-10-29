# mesh_generator.py
import gmsh
import config


def generate_mesh(show_gui=True):
    print("[1/5] Инициализация Gmsh...")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("SandwichPanel_3D")

    # Использование ядра OpenCASCADE для работы с 3D-солидами
    occ = gmsh.model.occ

    # 1. Создаем три сплошных слоя панели (блоки)
    # occ.addBox(x, y, z, dx, dy, dz, tag)
    layer1 = occ.addBox(0, 0, 0,
                        config.A, config.B, config.T_OUTER)

    layer2 = occ.addBox(0, 0, config.T_OUTER,
                        config.A, config.B, config.T_CORE)

    layer3 = occ.addBox(0, 0, config.T_OUTER + config.T_CORE,
                        config.A, config.B, config.T_INNER)

    # 2. Создаем объемы-инструменты для вырезания (окон и дверей)
    cutout_tools = []
    for (x, y, w, h) in config.CUTOUTS:
        # Вырез делаем чуть толще панели, чтобы гарантированно пробить насквозь
        z_start = -0.1
        z_depth = config.T_TOTAL + 0.2
        box = occ.addBox(x, y, z_start, w, h, z_depth)
        cutout_tools.append((3, box))  # 3 означает 3D-объект

    # Объединяем инструменты вырезов в один список
    tools = cutout_tools

    # 3. Вырезаем отверстия из каждого слоя отдельно,
    # чтобы слои остались независимыми объемами (для назначения разных материалов)
    print("[2/5] Булевы операции: вырезание окон...")

    # cut([(dim, tag)], tools, removeObject=True, removeTool=False)
    panel_out, _ = occ.cut([(3, layer1)], tools, removeTool=False)
    panel_core, _ = occ.cut([(3, layer2)], tools, removeTool=False)
    # На последнем слое удаляем инструменты-вырезы (removeTool=True)
    # ...
    panel_in, _ = occ.cut([(3, layer3)], tools, removeTool=True)

    occ.removeAllDuplicates()  # <--- ЭТА СТРОКА СКЛЕИВАЕТ СОПРИКАСАЮЩИЕСЯ ГРАНИ СЛОЕВ
    occ.synchronize()

    # 4. Назначаем физические группы
    # ...

    # 4. Назначаем физические группы (чтобы МКЭ решатель понимал, где какой материал)
    # panel_out[0][1] - возвращает тег созданного объема
    gmsh.model.addPhysicalGroup(3, [panel_out[0][1]], name="outer_skin")
    gmsh.model.addPhysicalGroup(3, [panel_core[0][1]], name="core")
    gmsh.model.addPhysicalGroup(3, [panel_in[0][1]], name="inner_skin")

    # TODO на следующих шагах: здесь же добавим физические группы для поверхностей
    # (жесткая заделка снизу и приложение нагрузки сверху).

    # 5. Генерация сетки
    print("[3/5] Генерация 3D сетки тетраэдров...")
    gmsh.option.setNumber("Mesh.MeshSizeMax", config.MESH_SIZE_MAX)
    gmsh.model.mesh.generate(3)

    # Сохраняем сетку в файл
    gmsh.write(config.MESH_FILE)
    print(f"[4/5] Сетка сохранена в файл {config.MESH_FILE}")

    # Запускаем графический интерфейс для проверки результата
    if show_gui:
        print("[5/5] Запуск интерфейса просмотра. Закрой окно Gmsh для продолжения.")
        gmsh.fltk.run()

    gmsh.finalize()


if __name__ == "__main__":
    generate_mesh(show_gui=True)