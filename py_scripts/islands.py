""" A Python script defining the extents of different island groups across the different PARTneR 2 nations (and some other pacific countries). 
Note apostropies in names are preserved. Capitalisation is not perserved. Spaces are replaced with underscores.

    Country keys are:
        cook_islands
        RMI_islands
        tonga
        vanuatu
        samoa
        tuvalu
        fiji (only Viti Levu)

"""

# UTM zones as EPSG CRS codes
crs_zone1S = 32701 # 180W to 174W
crs_zone2S = 32702 # 174W to 168W
crs_zone3S = 32703 # 168W to 162W
crs_zone4S = 32704 # 162W to 156W

crs_zone58S = 32758 # 162E to 168E
crs_zone59S = 32759 # 168E to 174E
crs_zone60S = 32760 # 174E to 180E

crs_zone57N = 32657 # 156E to 162E
crs_zone58N = 32658 # 162E to 168E
crs_zone59N = 32659 # 168E to 174E




## Function to return any of the country dictionaries by name
def get_island(country_name: str) ->dict:
    countries = {"cook_islands": cook_islands,
           "republic_of_marshall_islands": republic_of_marshall_islands,
           "samoa": samoa,
           "tonga": tonga,
           "tuvalu": tuvalu,
           "vanuatu": vanuatu,
                "fiji": fiji}
    country_function = countries[country_name]
    return country_function()


################################################################
# Cook Islands
def cook_islands() :
    cook_islands = {}

    cook_islands["mangaia"] = {"crs": crs_zone4S,
                         "lats": [-21.88, -21.88, -21.97, -21.97],
                         "lons": [-157.98, -157.85, -157.85, -157.98]}

    cook_islands["rarotonga"] = {"crs": crs_zone4S,
                      "lats": [-21.18, -21.18, -21.29, -21.29],
                        "lons": [-159.85, -159.70, -159.70, -159.85]}

    cook_islands["mauke"] = {"crs": crs_zone4S,
                        "lats": [-20.12, -20.12, -20.20, -20.20],
                        "lons": [-157.37, -157.31, -157.31, -157.37]}

    cook_islands["mitiaro"] = {"crs": crs_zone4S,
                            "lats": [-19.83, -19.83, -19.90, -19.90],
                            "lons": [-157.73, -157.63, -157.63, -157.73]}

    cook_islands["atiu"] = {"crs": crs_zone4S,
                        "lats": [-19.79, -19.79, -20.06, -20.06],
                        "lons": [-158.32, -158.06, -158.06, -158.32]}

    cook_islands["aitutaki_manuae_islands"] = {"crs": crs_zone4S,
                        "lats": [-18.81, -18.81, -19.30, -19.30],
                        "lons": [-159.86, -158.90, -158.90, -159.86]}

    cook_islands["palmerston"] = {"crs": crs_zone3S,
                        "lats": [-17.98, -17.98, -18.11, -18.11],
                        "lons": [-163.21, -163.10, -163.10, -163.21]}

    cook_islands["suwarrow"] = {"crs": crs_zone3S,
                        "lats": [-13.19, -13.19, -13.35, -13.35],
                        "lons": [-163.21, -163.03, -163.03, -163.21]}

    cook_islands["nassau"] = {"crs": crs_zone3S,
                        "lats": [-11.55, -11.55, -11.57, -11.57],
                        "lons": [-165.43, -165.40, -165.40, -165.43]}

    cook_islands["pukapuka"] = {"crs": crs_zone3S,
                        "lats": [-10.84, -10.84, -10.93, -10.93],
                        "lons": [-165.94, -165.82, -165.82, -165.94]}

    cook_islands["manihiki"] = {"crs": crs_zone4S,
                        "lats": [-10.35, -10.35, -10.49, -10.49],
                        "lons": [-161.06, -160.93, -160.93, -161.06]}

    cook_islands["rakahanga"] = {"crs": crs_zone4S,
                        "lats": [-9.99, -9.99, -10.05, -10.05],
                        "lons": [-161.11, -161.07, -161.07, -161.11]}

    cook_islands["penrhyn"] = {"crs": crs_zone4S,
                        "lats": [-8.90, -8.90, -9.10, -9.10],
                        "lons": [-158.07, -157.87, -157.87, -158.07]}
    
    return cook_islands




##############################################################
# Republic of Marshall Islands
def republic_of_marshall_islands() -> dict: 
    RMI_islands = {}

    ## Ralik Island chain
    RMI_islands["bikini"] = {"crs": crs_zone58N,
                         "lats": [11.75, 11.75, 11.46, 11.46],
                         "lons": [165.17, 165.61, 165.61, 165.17]}

    RMI_islands["rongelap_ailinginae_rongerik_atolls"] = {"crs": crs_zone58N,
                           "lats": [11.07, 11.07, 11.50, 11.50],
                           "lons": [166.25, 167.54, 167.54, 166.25]}

    RMI_islands["wotho"] = {"crs": crs_zone58N,
                        "lats": [10.19, 10.19, 10.00, 10.00],
                        "lons": [165.91, 166.05, 166.05, 165.91]}

    RMI_islands["kwajalein"] = {"crs": crs_zone58N,
                            "lats": [9.43, 9.43, 8.65, 8.65],
                            "lons": [166.75, 167.85, 167.85, 166.75]}

    RMI_islands["ujae"] = {"crs": crs_zone58N,
                        "lats": [9.24, 9.24, 8.88, 8.88],
                        "lons": [165.49, 165.82, 165.82, 165.49]}

    RMI_islands["lae"] = {"crs": crs_zone58N,
                      "lats": [8.98, 8.98, 8.91, 8.91],
                      "lons": [166.28, 166.20, 166.20, 166.28]}

    RMI_islands["lib"] = {"crs": crs_zone58N,
                      "lats": [8.33, 8.33, 8.30, 8.30],
                      "lons": [167.36, 167.39, 167.39, 167.36]}

    RMI_islands["namu"] = {"crs": crs_zone59N,
                       "lats": [8.21, 8.21, 7.73, 7.73],
                       "lons": [167.96, 168.32, 168.32, 167.96]}

    RMI_islands["jabwot"] = {"crs": crs_zone59N,
                        "lats": [7.76, 7.76, 7.74, 7.74],
                        "lons": [168.97, 168.99, 168.99, 168.97]}

    RMI_islands["ailinglaplap"] = {"crs": crs_zone59N,
                               "lats": [7.62, 7.62, 7.24, 7.24],
                               "lons": [168.52, 168.99, 168.99, 168.52]}

    RMI_islands["jaluit"] = {"crs": crs_zone59N,
                         "lats": [6.62, 6.62, 5.77, 5.77],
                         "lons": [169.36, 169.74, 169.74, 169.36]}

    RMI_islands["kili"] = {"crs": crs_zone59N,
                       "lats": [5.66, 5.66, 5.62, 5.62],
                       "lons": [169.10, 169.15, 169.15, 169.10]}

    RMI_islands["ebon"] = {"crs": crs_zone59N,
                       "lats": [4.69, 4.69, 4.56, 4.56],
                       "lons": [168.63, 168.79, 168.79, 168.63]}


    ## Ratak Island chain
    RMI_islands["utirik_taka_atolls"] = {"crs": crs_zone59N,
                         "lats": [11.34, 11.34, 11.06, 11.06],
                         "lons": [169.50, 169.89, 169.89, 169.50]}

    RMI_islands["ailuk"] = {"crs": crs_zone59N,
                        "lats": [10.48, 10.48, 10.20, 10.20],
                        "lons": [169.83, 170.00, 170.00, 169.83]}

    RMI_islands["mejit"] = {"crs": crs_zone59N,
                        "lats": [10.31, 10.31, 10.26, 10.26],
                        "lons": [170.85, 170.88, 170.88, 170.85]}

    RMI_islands["likiep_jemo_atolls "] = {"crs": crs_zone59N,
                         "lats": [10.12, 10.12, 9.78, 9.78],
                         "lons": [168.97, 169.58, 169.58, 168.97]}

    RMI_islands["wotje_erikub_atolls"] = {"crs": crs_zone59N,
                        "lats": [9.57, 9.57, 8.99, 8.99],
                        "lons": [169.78, 170.29, 170.29, 169.78]}

    RMI_islands["maloelap"] = {"crs": crs_zone59N,
                           "lats": [8.92, 8.92, 8.47, 8.47],
                           "lons": [170.81, 171.27, 171.27, 170.81]}

    RMI_islands["aur"] = {"crs": crs_zone59N,
                      "lats": [8.37, 8.37, 8.12, 8.12],
                      "lons": [171.00, 171.19, 171.19, 171.00]}

    RMI_islands["majuro"] = {"crs": crs_zone59N,
                         "lats": [7.23, 7.23, 7.03, 7.03],
                         "lons": [171, 171.43, 171.43, 171]}

    RMI_islands["arno"] = {"crs": crs_zone59N,
                       "lats": [7.31, 7.31, 6.95, 6.95],
                       "lons": [171.52, 171.97, 171.97, 171.52]}

    RMI_islands["mili_knox_atolls"] = {"crs": crs_zone59N,
                       "lats": [6.28, 6.28, 5.86, 5.86],
                       "lons": [171.69, 172.19, 172.19, 171.69]}

    ## Outlier islands
    RMI_islands["namorik"] = {"crs": crs_zone59N,
                       "lats": [5.65, 5.65, 5.58, 5.58],
                       "lons": [168.07, 168.14, 168.14, 168.07]}

    RMI_islands["bikar"] = {"crs": crs_zone59N,
                       "lats": [12.29, 12.29, 12.18, 12.18],
                       "lons": [170.06, 170.15, 170.15, 170.06]}

    RMI_islands["bokak"] = {"crs": crs_zone59N,
                         "lats": [14.73, 14.73, 14.54, 14.54],
                         "lons": [168.88, 169.03, 169.03, 168.88]}

    RMI_islands["enewetak"] = {"crs": crs_zone58N,
                           "lats": [11.68, 11.68, 11.31, 11.31],
                           "lons": [162.02, 162.45, 162.45, 162.02]}

    RMI_islands["ujelang"] = {"crs": crs_zone57N,
                         "lats": [9.88, 9.88, 9.76, 9.76],
                         "lons": [160.79, 160.99, 160.99, 160.79]}
    
    return RMI_islands


######################################################################
# Tonga
def tonga() -> dict: 
    tonga = {}


    tonga["minerva"] = {"crs": crs_zone1S, 
                         "lats": [-23.59, -23.59, -23.98, -23.98],
                         "lons": [-179.19, -178.86, -178.86, -179.19]}

    tonga["'ata"] = {"crs": crs_zone1S,
                      "lats": [-22.32, -22.32, -22.36, -22.36],
                        "lons": [-176.23, -176.18, -176.18, -176.23]}

    tonga["'eua"] = {"crs": crs_zone1S,
                        "lats": [-21.26, -21.26, -21.49, -21.49],
                        "lons": [-174.99, -174.88, -174.88, -174.99]}

    tonga["tongatapu"] = {"crs": crs_zone1S,
                            "lats": [-20.98, -20.98, -21.29, -21.29],
                            "lons": [-174.95, -175.37, -175.37, -174.95]}

    tonga["hunga_tonga-hunga_ha'apai_islands"] = {"crs": crs_zone1S,
                        "lats": [-20.52, -20.52, -20.58, -20.58],
                        "lons": [-175.43, -175.34, -175.34, -175.43]}

    tonga["ha'apai_group"] = {"crs": crs_zone1S,
                        "lats": [-19.53, -19.53, -20.64, -20.64],
                        "lons": [-175.12, -174.21, -174.21, -175.12]}

    tonga["vava'u_group"] = {"crs": crs_zone1S,
                        "lats": [-18.01, -18.01, -18.99, -18.99],
                        "lons": [-174.86, -173.89, -173.89, -174.86]}

    tonga["niuatoputapu_tafahi_islands"] = {"crs": crs_zone2S,
                        "lats": [-15.82, -15.82, -16.00, -16.00],
                        "lons": [-173.85, -173.71, -173.71, -173.85]}

    tonga["niuafo'ou"] = {"crs": crs_zone1S,
                        "lats": [-15.55, -15.55, -15.65, -15.65],
                        "lons": [-175.68, -175.59, -175.59, -175.68]}
    
    return tonga



#######################################################################
# Vanuatu
def vanuatu() -> dict:
    vanuatu = {}

    ## Espiritu
    vanuatu["malampa_province"] = {"crs": crs_zone58S,
                            "lats": [-15.83, -16.79, -16.56, -16.56, -16.07, -16.07],
                            "lons": [167.12, 167.12, 167.99, 168.42, 168.42, 167.99]}

    vanuatu["penama_province"] = {"crs": crs_zone58S,
                        "lats": [-14.89, -14.89, -16.03, -16.03],
                        "lons": [167.66, 168.29, 168.29, 167.66]}

    vanuatu["tafea_province"] = {"crs": crs_zone59S,
                        "lats": [-18.61, -18.61, -20.28, -20.28],
                        "lons": [168.96, 170.25, 170.25, 168.96]}

    vanuatu["torba_province "] = {"crs": crs_zone58S,
                        "lats": [-13.05, -13.05, -14.51, -14.51],
                        "lons": [166.50, 168.09, 168.09, 166.50]}

    vanuatu["shefa_province "] = {"crs": crs_zone59S,
                        "lats": [-16.55, -16.55, -17.84, -17.84],
                        "lons": [168.10, 168.68, 168.68, 168.10]}

    vanuatu["sanma_province"] = {"crs": crs_zone58S,
                         "lats": [-14.63, -14.63, -15.77, -15.77],
                         "lons": [166.53, 167.31, 167.31, 166.53]}
    
    return vanuatu


#####################################################################
# Samoa
def samoa() -> dict:
    samoa = {}

    samoa["savai'i"] = {"crs": crs_zone2S,
                          "lats": [-13.4, -13.83, -13.83, -13.4],
                          "lons": [-172.82, -172.82, -172.15, -172.15]}

    samoa["upolu"] = {"crs": crs_zone2S,
                        "lats": [-13.75, -14.08, -14.08, -13.75],
                        "lons": [-172.14, -172.14, -171.35, -171.35]}
    
    return samoa

#####################################################################
# Tuvalu
def tuvalu() -> dict: 
    tuvalu = {}

    tuvalu["vaitupu"] = {"crs": crs_zone60S,
                          "lats": [-7.44,-7.5,  -7.5, -7.44],
                          "lons": [178.6, 178.6, 178.71, 178.71]}

    tuvalu["nukufetau"] = {"crs": crs_zone60S,
                          "lats": [-9.3, -9.45, -9.45, -9.3],
                          "lons": [179.77, 179.77, 179.9, 179.9]}

    tuvalu["nanumanga"] = {"crs": crs_zone60S,
                          "lats": [-5.6, -6.31, -6.31, -5.6],
                          "lons": [176.3, 176.3, 176.35, 176.35]}

    tuvalu["nanumea"] = {"crs": crs_zone60S,
                          "lats": [-5.6, -5.72, -5.72, -5.6],
                          "lons": [176, 176, 176.2, 176.2]}

    tuvalu["niulakita"] = {"crs": crs_zone60S,
                          "lats": [-10.78, -10.8, -10.8, -10.78],
                          "lons": [179.45, 179.45, 179.5, 179.5]}

    tuvalu["niutao"] = {"crs": crs_zone60S,
                          "lats": [-6.08, -6.13, -6.13, -6.08],
                          "lons": [177.3, 177.3, 177.4, 177.4]}

    tuvalu["nui"] = {"crs": crs_zone60S,
                     "lats": [-7.1, -7.3, -7.3, -7.1],
                     "lons": [177.1, 177.1, 177.2, 177.2]}

    tuvalu["nukulaelae"] = {"crs": crs_zone60S,
                           "lats": [-7.9, -8.1, -8.1, -7.9],
                           "lons": [178.2, 178.2, 178.5, 178.5]}

    tuvalu["funafuti"] = {"crs": crs_zone60S,
                          "lats": [-8.4, -8.7, -8.7, -8.4],
                          "lons": [179, 179, 179.3, 179.3]}
    
    return tuvalu



##############################################################
# Republic of Marshall Islands
def fiji() -> dict:
    fiji = {}


    fiji["viti-levu"] = {"crs": crs_zone60S,
                         "lats": [-17.2, -17.2, -18.7, -18.7],
                         "lons": [175, 179, 179, 175]}
    return fiji
