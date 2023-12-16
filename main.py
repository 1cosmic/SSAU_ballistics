from docxtpl import DocxTemplate
from calculates import *

document = DocxTemplate("test.docx")

# Consts.
g_const = 9.801

# INPUT YOU VALUES ON FERE:
rocket = {
    'name': "Тор-Эджена Д",

    # First part of Rocket:
    'm': (48.6, 6.92),
    'mt': (45, 6.15),
    'mk': (3.6, 0.77),
    'P0': (765, "-"),
    'P_vacuum': (860, 71),
    'P_ud': (250 * g_const, "-"),
    'P_ud_vacuum': (284 * g_const, 320 * g_const),
    'D_midel': (3.5, 3.5),


    # Other parameters:
    'm_pn': 1,
    'type_of_fuel': "Аэрозин + тетраоксид",
}

ballistic_data = {
    'H': 210,
    'i': 51.6,
    'fi_0': 31.2,
}

ballistic_pars = ballistic_parameters(rocket)
for key in ballistic_pars:
    print(key, ' ', ballistic_pars[key])

#
# document.render(content)
# document.save("saved_test.docx")