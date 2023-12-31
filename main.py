import docx
from docxtpl import DocxTemplate
from calculates import *

document = DocxTemplate("test.docx")

# Consts.
g_const = 9.801

# INPUT YOU VALUES ON FERE:
# rocket = {
#     'name': "Тор-Эджена Д",
#
#     # First part of Rocket:
#     'm': (48.6, 6.92),
#     'mt': (45, 6.15),
#     'mk': (3.6, 0.77),
#     'P0': (765, "-"),
#     'P_vacuum': (860, 71),
#     'P_ud': (250 * g_const, "-"),
#     'P_ud_vacuum': (284 * g_const, 320 * g_const),
#     'D_midel': (3.5, 3.5),
#
#
#     # Other parameters:
#     'm_pn': 1,
#     'type_of_fuel': "Аэрозин + тетраоксид",
# }
#
# ballistic_data = {
#     'H': 210,
#     'i': 51.6,
#     'fi_0': 31.2,
# }

# Цыцаров
m1 = 95
m2 = 26
mt1 = 88
mt2 = 24
mk1 = 7
mk2 = 2
P0_1 = 2100
P0_2 = "-"
P_vacuum_1 = 2415
P_vacuum_2 = 900
P_ud_1 = 2746
P_ud_2 = "-"
P_ud_vacuum_1 = 3158
P_ud_vacuum_2 = 2835
D_midel_1 = 3
D_midel_2 = 1.9

m_pn = 1
rocket_name = "Вега"
user_name = "Попов Л.В."
number_exs = "9"

ballistic_data = {
    'H': 200,
    'i': 56,
    'i_target': 51.6,
    'fi_0': 51.9,
    'target_H': 415,
}

# # Тест
# m1 = 95
# m2 = 26
# mt1 = 88
# mt2 = 24
# mk1 = 7
# mk2 = 2
# P0_1 = 2100
# P0_2 = "-"
# P_vacuum_1 = 2415
# P_vacuum_2 = 900
# P_ud_1 = 2746
# P_ud_2 = "-"
# P_ud_vacuum_1 = 3158
# P_ud_vacuum_2 = 2835
# D_midel_1 = 3
# D_midel_2 = 1.9
# m_pn = 1
#
# ballistic_data = {
#     'H': 210,
#     'i': 20,
#     'i_target': 51.6,
#     'fi_0': 13.9,
#     'target_H': 415,
# }

rocket = {
    'name': rocket_name,
    'user_name': user_name,
    'number_exs': number_exs,

    # First part of Rocket:
    'm': (m1, m2),
    'mt': (mt1, mt2),
    'mk': (mk1, mk2),
    'P0': (P0_1, P0_2),
    'P_vacuum': (P_vacuum_1, P_vacuum_2),
    'P_ud': (P_ud_1, P_ud_2),
    'P_ud_vacuum': (P_ud_vacuum_1, P_ud_vacuum_2),
    'D_midel': (D_midel_1, D_midel_2),

    'm1': m1,
    'm2': m2,
    'mt1': mt1,
    'mt2': mt2,
    'mk1': mk1,
    'mk2': mk2,
    'P0_1': P0_1,
    'P0_2': P0_2,
    'P_vacuum_1': P_vacuum_1,
    'P_vacuum_2': P_vacuum_2,
    'P_ud_1': P_ud_1,
    'P_ud_2': P_ud_2,
    'P_ud_vacuum_1': P_ud_vacuum_1,
    'P_ud_vacuum_2': P_ud_vacuum_2,
    'D_midel_1': D_midel_1,
    'D_midel_2': D_midel_2,

    # Other parameters:
    'm_pn': m_pn,
    'type_of_fuel': "Аэрозин + тетраоксид",
}


data_ballistic = ballistic_parameters(rocket)

with open("data_ballistic.txt", "w") as file:
    for key in data_ballistic:
        text = f"{key}: {data_ballistic[key]}"
        file.write(text + '\n')
        # print(text)

data_Vx = calcs_Vx(ballistic_data)
data_corr_mass = correct_mass_fuel()
data_iner = swap_start_to_iner()

all_data = data_ballistic | data_corr_mass | data_Vx | data_iner

clean_data_for_render = {}
for key in all_data:
    if isinstance(all_data[key], tuple) != True:
        if isinstance(all_data[key], float) and key != 'tetta0':
            all_data[key] = (str(round(all_data[key], 3))).replace('.', ',')


    # print(key, all_data[key])

document.render(all_data)
document.save(rocket['user_name'] + ".docx")