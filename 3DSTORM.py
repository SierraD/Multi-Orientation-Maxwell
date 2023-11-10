import numpy as np
import pandas as pd

Name1 = input("XY Data file name: ")
Name1=str(Name1)
Name2 = input("XZ Data file name: ")
Name2=str(Name2)
XYPIX = input("XY pixel size [nm]: ")
XYPIX=int(XYPIX)
ZPIX = input("Z pixel size [nm]: ")
ZPIX=int(ZPIX)
XYSEARCH = input("XY search radius [nm]: ")
XYSEARCH=int(XYSEARCH)
ZSEARCH = input("Z search radius [nm]: ")
ZSEARCH=int(ZSEARCH)
print("Please wait...")

def Set_Units(XY, Z, name_xy, name_xz):
    """Open and read ThunderSTORM Analysis files and save the data with standard formatting.
    
    Args: 
    XY: The in-plane pixel size which comprises the 3D voxel, in nanometers.
    Z: The in-depth pixel length which comprises the 3D voxel, in nanometers. 
    name_xy: The name of the ThunderSTORM .csv file which contains the data obtained from the XY Orientation.
    name_xz: The name of the ThunderSTORM .csv file which contains the data obtained from the XZ Orientation.
    
    Returns:
    Data: A dataset containing all of the ThunderSTORM values, separated into standart formatted columns.
    """
    df_xy = pd.read_csv(name_xy)
    df_xz = pd.read_csv(name_xz)
    Data = {"X_XY":df_xy["x [nm]"]*XY, 
            "Y_XY":df_xy["y [nm]"]*XY, 
            "Z_XY":df_xy["frame"]*Z, 
            "X_XZ":df_xz["x [nm]"]*XY, 
            "Y_XZ":df_xz["frame"]*XY, 
            "Z_XZ":df_xz["y [nm]"]*Z, 
            "Uncertainty_XY":df_xy["uncertainty [nm]"]*XY, 
            "UncertaintyX_XZ":df_xz["uncertainty [nm]"]*XY, 
            "UncertaintyZ_XZ":df_xz["uncertainty [nm]"]*Z,
            "Intensity_XY":df_xy["intensity [photon]"], 
            "Intensity_XZ":df_xz["intensity [photon]"], 
            "Sigma_XY":df_xy["sigma [nm]"]*XY, 
            "SigmaX_XZ":df_xz["sigma [nm]"]*XY, 
            "SigmaZ_XZ":df_xz["sigma [nm]"]*Z}
    return Data

def Find_Overlap(PM_y, PM_z, data):
    """For each XZ localization, search all XY localizations for a value which overlaps within the 3D Uncertainty. 
    
    Args:
    data: A dictionary of all X, Y, Z localizations and associated uncertainties obtained from the XY and XZ data.
    PM_y: The integer value of the uncertainty in Y, which is the step size for the XZ data.
    PM_z: The integer value of the uncertainty in Z, which is the step size for the XY data.
    
    Returns:
    overlapped_data: A dictionary of particle data where overlap between the XY and XZ data occurs in standard formatting.
    """
    XY_Index = []
    XZ_Index = []
    for i in range(0, len(data["X_XZ"])):
        count = 0
        PM_x = data["UncertaintyX_XZ"][i]
        for j in range(0, len(data["X_XY"])):
            if (data["X_XZ"][i] - PM_x < data["X_XY"][j] < data["X_XZ"][i] + PM_x):
                if (data["Y_XZ"][i] - PM_y < data["Y_XY"][j] < data["Y_XZ"][i] + PM_y):
                    if (data["Z_XZ"][i] - PM_z < data["Z_XY"][j] < data["Z_XZ"][i] + PM_z):
                        count = count+1
            if (count>0):
                XY_Index.append(j)
                XZ_Index.append(i)
                count=0           
    overlap = [[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    for k in XY_Index:
        overlap[0].append(data["X_XY"][k]) 
        overlap[1].append(data["Y_XY"][k]) 
        overlap[2].append(data["Z_XY"][k])
        overlap[3].append(data["Uncertainty_XY"][k]) 
        overlap[4].append(data["Intensity_XY"][k])
        overlap[5].append(data["Sigma_XY"][k]) 
    for l in XZ_Index:
        overlap[6].append(data["X_XZ"][l]) 
        overlap[7].append(data["Y_XZ"][l]) 
        overlap[8].append(data["Z_XZ"][l])
        overlap[9].append(data["UncertaintyX_XZ"][l]) 
        overlap[10].append(data["UncertaintyZ_XZ"][l]) 
        overlap[11].append(data["Intensity_XZ"][l]) 
        overlap[12].append(data["SigmaX_XZ"][l]) 
        overlap[13].append(data["SigmaZ_XZ"][l]) 
        overlapped_data = {"X_XY":overlap[0], "Y_XY":overlap[1], "Z_XY":overlap[2],
                 "Uncertainty_XY":overlap[3], 
                 "X_XZ":overlap[6], "Y_XZ":overlap[7], "Z_XZ":overlap[8], 
                 "UncertaintyX_XZ":overlap[9], "UncertaintyZ_XZ":overlap[10], 
                 "Intensity_XY":overlap[4], "Intensity_XZ":overlap[11], 
                 "Sigma_XY":overlap[5], "SigmaX_XZ":overlap[12], "SigmaZ_XZ":overlap[13]}
    return overlapped_data

def Uncertainty_Calculation(data):
    """Using the data acquired from the overlap calculation, which includes can include many duplicate points, filter through all of 
    the duplicates and only save the overlapped values with the smallest uncertainty.
    
    Args: 
    data: A dictionary of particle data where overlap between the XY and XZ data occurs, in standard formatting.
    
    Returns:
    data_updated: A dictionary which includes the overlapped data without any duplicate points, as selected by uncertainty, in 
    standard formatting.
    """
    for s in range(0, 2):
        if s == 0:
            orientation = "XZ"
        else:
            orientation = "XY"
            data = data_updated
        DF = pd.DataFrame(data=data)
        Groupby = DF.groupby("X_"+str(orientation)).count().reset_index()
        all_values = Groupby["X_"+str(orientation)]
        if orientation == "XY":
            other_orientation = "XZ"
            uncert_orientation = "UncertaintyZ_XZ"
        elif orientation == "XZ":
            other_orientation = "XY"
            uncert_orientation = "Uncertainty_XY"
        certain_points = [[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
        for i in range(0, len(np.unique(data["X_"+str(orientation)]))):
            locations = np.where(data["X_"+str(orientation)]==all_values[i])[0] 
            values_X = []
            values_Y = []
            values_Z = []
            values_Uncert = []
            for j in range(0, len(locations)):
                values_X.append(data["X_"+str(other_orientation)][locations[j]])
                values_Y.append(data["Y_"+str(other_orientation)][locations[j]])
                values_Z.append(data["Z_"+str(other_orientation)][locations[j]])
                values_Uncert.append(data[str(uncert_orientation)][locations[j]])
            min_uncert = np.where(values_Uncert == min(values_Uncert))
            index_uncert = locations[min_uncert[0]]              
            certain_points[0].append(data["X_XY"][index_uncert[0]])
            certain_points[1].append(data["Y_XY"][index_uncert[0]])
            certain_points[2].append(data["Z_XY"][index_uncert[0]])
            certain_points[3].append(data["X_XZ"][index_uncert[0]])
            certain_points[4].append(data["Y_XZ"][index_uncert[0]])
            certain_points[5].append(data["Z_XZ"][index_uncert[0]])
            certain_points[6].append(data["Uncertainty_XY"][index_uncert[0]])
            certain_points[7].append(data["UncertaintyX_XZ"][index_uncert[0]])
            certain_points[8].append(data["UncertaintyZ_XZ"][index_uncert[0]])
            certain_points[9].append(data["Intensity_XY"][index_uncert[0]])
            certain_points[10].append(data["Intensity_XZ"][index_uncert[0]])
            certain_points[11].append(data["Sigma_XY"][index_uncert[0]])
            certain_points[12].append(data["SigmaX_XZ"][index_uncert[0]])
            certain_points[13].append(data["SigmaZ_XZ"][index_uncert[0]])
            data_updated = {"X_XY":certain_points[0], "Y_XY":certain_points[1], "Z_XY":certain_points[2], 
                          "X_XZ":certain_points[3], "Y_XZ":certain_points[4], "Z_XZ":certain_points[5], 
                          "Uncertainty_XY":certain_points[6], 
                          "UncertaintyX_XZ":certain_points[7], "UncertaintyZ_XZ":certain_points[8],
                          "Intensity_XY":certain_points[9], "Intensity_XZ":certain_points[10], 
                          "Sigma_XY":certain_points[11], "SigmaX_XZ":certain_points[12], "SigmaZ_XZ":certain_points[12]}
    return data_updated

def Proximity_Calculation(data):
    """Using the data acquired from the overlap calculation, which includes can include many duplicate points, filter through all 
    of the duplicates and only save the overlapped values with the closest 3D proximity.
    
    Args: 
    data: A dictionary of particle data where overlap between the XY and XZ data occurs, in standard formatting.
    
    Returns:
    data_updated: A dictionary which includes the overlapped data without any duplicate points, as selected by proximity, in 
    standard formatting.
    """
    for s in range(0, 2):
        if s == 0:
            orientation = "XZ"
        else:
            orientation = "XY"
            data = data_updated
        DF = pd.DataFrame(data=data)
        Groupby = DF.groupby("X_"+str(orientation)).count().reset_index()
        all_values = Groupby["X_"+str(orientation)]
        if orientation == "XY":
            other_orientation = "XZ"
            uncert_orientation = "UncertaintyZ_XZ"
        elif orientation == "XZ":
            other_orientation = "XY"
            uncert_orientation = "Uncertainty_XY"
        close_points = [[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
        for i in range(0, len(np.unique(data["X_"+str(orientation)]))):
            locations = np.where(data["X_"+str(orientation)]==all_values[i])[0]
            values_X = []
            values_Y = []
            values_Z = []
            for j in range(0, len(locations)):
                values_X.append(data["X_"+str(other_orientation)][locations[j]])
                values_Y.append(data["Y_"+str(other_orientation)][locations[j]])
                values_Z.append(data["Z_"+str(other_orientation)][locations[j]])
            x_dist = []
            y_dist = []
            z_dist = []
            for k in range(0, len(locations)):
                x_dist.append(np.abs(data["X_"+str(orientation)][locations[0]] - values_X[k]))
                y_dist.append(np.abs(data["Y_"+str(orientation)][locations[0]] - values_Y[k]))
                z_dist.append(np.abs(data["X_"+str(orientation)][locations[0]] - values_Z[k]))
            pythagoras = []
            for l in range(0, len(locations)):
                pythagoras.append(np.sqrt(x_dist[l]**2 + y_dist[l]**2 + z_dist[l]**2))     
            min_close = np.where(pythagoras == min(pythagoras))
            index_close = locations[min_close[0]]
            close_points[0].append(data["X_XY"][index_close[0]])
            close_points[1].append(data["Y_XY"][index_close[0]])
            close_points[2].append(data["Z_XY"][index_close[0]])
            close_points[3].append(data["X_XZ"][index_close[0]])
            close_points[4].append(data["Y_XZ"][index_close[0]])
            close_points[5].append(data["Z_XZ"][index_close[0]])
            close_points[6].append(data["Uncertainty_XY"][index_close[0]])
            close_points[7].append(data["UncertaintyX_XZ"][index_close[0]])
            close_points[8].append(data["UncertaintyZ_XZ"][index_close[0]])
            close_points[9].append(data["Intensity_XY"][index_close[0]])
            close_points[10].append(data["Intensity_XZ"][index_close[0]])
            close_points[11].append(data["Sigma_XY"][index_close[0]])
            close_points[12].append(data["SigmaX_XZ"][index_close[0]])
            close_points[13].append(data["SigmaZ_XZ"][index_close[0]])
            data_updated = {"X_XY":close_points[0], "Y_XY":close_points[1], "Z_XY":close_points[2], 
                          "X_XZ":close_points[3], "Y_XZ":close_points[4], "Z_XZ":close_points[5], 
                          "Uncertainty_XY":close_points[6], 
                          "UncertaintyX_XZ":close_points[7], "UncertaintyZ_XZ":close_points[8],
                          "Intensity_XY":close_points[9], "Intensity_XZ":close_points[10], 
                          "Sigma_XY":close_points[11], "SigmaX_XZ":close_points[12], "SigmaZ_XZ":close_points[13]}
    return data_updated

def Isolated_Points(data):
    """Using the data acquired from either the uncertainty or the proximity calculation, filter through all of the localizations 
    and determine if any of them are within 3D uncertainty of others, indicating that it is multiple measurements from the same 
    particle. If a group of multiple measurements is found, select and save only the measurement with the smallest uncertainty.
    
    Args:
    data: A dictionary of XY and XZ overlapped points, containing no duplicates after calculation from the uncertainty method 
    or the proximity method, in standard formatting.
    
    Returns:
    data_final: A dictionary of data which includes the overlapped data without any duplicate points or multiple measurements of 
    the same localization, in standard formatting.
    """
    Isolated = []
    Grouped = []
    for i in range(0, len(data["X_XY"])):
        X = data["X_XY"][i]
        Y = data["Y_XY"][i]
        Z = data["Z_XZ"][i]
        L_X = X - data["Uncertainty_XY"][i]
        R_X = X + data["Uncertainty_XY"][i]
        L_Y = Y - data["Uncertainty_XY"][i]
        R_Y = Y + data["Uncertainty_XY"][i]
        L_Z = Z - data["UncertaintyZ_XZ"][i]
        R_Z = Z + data["UncertaintyZ_XZ"][i]
        indexes = []
        for j in range(0, len(data["X_XY"])):
            SX = data["X_XY"][j]
            SY = data["Y_XY"][j]
            SZ = data["Z_XZ"][j]
            SL_X = SX - data["Uncertainty_XY"][j]
            SR_X = SX + data["Uncertainty_XY"][j]
            SL_Y = SY - data["Uncertainty_XY"][j]
            SR_Y = SY + data["Uncertainty_XY"][j]
            SL_Z = SZ - data["UncertaintyZ_XZ"][j]
            SR_Z = SZ + data["UncertaintyZ_XZ"][j]
            if ((L_X<SL_X<R_X)or(L_X<SR_X<R_X))or((SL_X>L_X)and(SR_X<R_X)):
                if ((L_Y<SL_Y<R_Y)or(L_Y<SR_Y<R_Y))or((SL_Y>L_Y)and(SR_Y<R_Y)):
                    if ((L_Z<SL_Z<R_Z)or(L_Z<SR_Z<R_Z))or((SL_Z>L_Z)and(SR_Z<R_Z)):
                        indexes.append(j)
        indexes.append(i)
        indexes.sort()
        if (len(indexes)==1):
            Isolated.append(indexes)
        else:
            Grouped.append(indexes)    
    Isolated = np.asarray(Isolated).flatten()
    Groups = []
    for item in Grouped:
        if item not in Groups:
            Groups.append(item)
    Grouped_Selection = []
    for k in range(0, len(Grouped)):
        Uncertainty_Number = []
        for l in range(0, len(Grouped[k])):
            Uncertainty_Number.append(data["Uncertainty_XY"][Grouped[k][l]]+data["UncertaintyZ_XZ"][Grouped[k][l]])
        Selection = np.where(Uncertainty_Number == min(Uncertainty_Number))
        Grouped_Selection.append(Grouped[k][Selection[0][0]])
    Particle_Indexes = np.sort(np.concatenate([Isolated, Grouped_Selection]))
    points = [[],[],[],[],[],[],[],[],[]]
    for m in Particle_Indexes:
        points[0].append(data["X_XY"][m])
        points[1].append(data["Y_XY"][m])
        points[2].append(data["Z_XZ"][m])
        points[3].append(data["Uncertainty_XY"][m])
        points[4].append(data["UncertaintyZ_XZ"][m])
        points[5].append(data["Intensity_XY"][m])
        points[6].append(data["Intensity_XZ"][m])
        points[7].append(data["Sigma_XY"][m])
        points[8].append(data["SigmaZ_XZ"][m])
    data_final = {"X [nm]":points[0], "Y [nm]":points[1], "Z [nm]":points[2], 
            "XY Uncertainty [nm]":points[3], "Z Uncertainty [nm]":points[4], 
            "XY Intensity [A.U.]":points[5], "Z Intensity [A.U.]":points[6], 
            "XY Sigma [nm]":points[7], "Z Sigma [nm]":points[8]}
    return data_final
Raw_Data = Set_Units(XYPIX, ZPIX, Name1, Name2)
Overlap_Data = Find_Overlap(XYSEARCH, ZSEARCH, Raw_Data)

Certain_Data = Uncertainty_Calculation(Overlap_Data)
Close_Data = Proximity_Calculation(Overlap_Data)
Final_Uncertainty = Isolated_Points(Certain_Data)
Final_Proximity = Isolated_Points(Close_Data)
Final_Uncertainty_File = pd.DataFrame(data=Final_Uncertainty)
Final_Uncertainty_File.to_csv("3DSTORM_UncertaintyMethod.csv", sep=',', index=False, encoding='utf-8')
Final_Proximity_File = pd.DataFrame(data=Final_Proximity)
Final_Proximity_File.to_csv("3DSTORM_ProximityMethod.csv", sep=',', index=False, encoding='utf-8')
print("Download Complete")