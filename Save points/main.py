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

m1 = 118.8
m2 = 28.4
mt1 = 113
mt2 = 26.4
mk1 = 5.8
mk2 = 2
P0_1 = 1915
P0_2 = "-"
P_vacuum_1 = 2078
P_vacuum_2 = 445
P_ud_1 = 270 * g_const
P_ud_2 = "-"
P_ud_vacuum_1 = 293 * g_const
P_ud_vacuum_2 = 316 * g_const
D_midel_1 = 4.2
D_midel_2 = 4

rocket = {
    'name': "Титан-2",
    'user_name': "Захаров И.М.",
    'number_exs': 3,

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
    'm_pn': 4,
    'type_of_fuel': "Аэрозин + тетраоксид",
}

ballistic_data = {
    'H': 190,
    'i': 63,
    'i_target': 51.6,
    'fi_0': 62.8,
    'target_H': 415,
}

data_ballistic = ballistic_parameters(rocket)
data_Vx = calcs_Vx(ballistic_data)
data_corr_mass = correct_mass_fuel()
data_iner = swap_start_to_iner()

all_data = data_ballistic | data_corr_mass | data_Vx | data_iner

#
# for key in all_data:
#     if isinstance(all_data[key], float) and key != 'tetta0':
#         all_data[key] = round(all_data[key], 3)
#
#     print(key, all_data[key])

document.render(all_data)
document.save("saved_test.docx")