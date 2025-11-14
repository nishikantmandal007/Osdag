"""

@Author:    Rutvik Joshi - Osdag Team, IIT Bombay [(P) rutvikjoshi63@gmail.com / 30005086@iitb.ac.in]
12.03.2025
Revised Design for GUI: Parth Karia - Osdag Team, IIT Bombay [30006096@iitb.ac.in]

@Module - Beam Design- Simply Supported member
           - Laterally Supported Beam [Moment + Shear]
           - Laterally Unsupported Beam [Moment + Shear]


@Reference(s): 1) IS 800: 2007, General construction in steel - Code of practice (Third revision)
               2) IS 808: 1989, Dimensions for hot rolled steel beam, column, channel, and angle sections and
                                it's subsequent revision(s)
               3) Design of Steel Structures by N. Subramanian (Fifth impression, 2019, Chapter 15)
               4) Limit State Design of Steel Structures by S K Duggal (second edition, Chapter 11)

other          8)
references     9)

"""
import logging
import math
import numpy as np
# Use the standard pyswarm PSO library
from pyswarm import pso

# --- Osdag Core Imports ---
from ...Common import *
from ...utils.common.material import *
from ...utils.common.load import Load
from ...utils.common.component import ISection, Material, Plate
from ...utils.common import is800_2007
from ...utils.common.common_calculation import *
from ...utils.common.Unsymmetrical_Section_Properties import Unsymmetrical_I_Section_Properties
from ..member import Member
from ...Report_functions import *
from osdag.cad.items.plate import Plate

# --- PyQt5 Imports (needed for base classes) ---
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QDialog
from PyQt5.QtCore import Qt

# --- Refactored Module Imports ---
from .custom_widgets import RangeInputDialog, My_ListWidget, My_ListWidgetItem, PopupDialog
from .pso_solver import Section
from .design_checks import PlateGirderLogic


# Constants for plate girder design
KEY_OVERALL_DEPTH_PG_CST = "Overall Depth (D) (mm)"


class PlateGirderWelded(Member, PlateGirderLogic):
    """
    Main class for Welded Plate Girder design.
    
    This class acts as the "controller" interfacing with the Osdag framework.
    It inherits core engineering logic from PlateGirderLogic.
    It uses GUI components from custom_widgets.
    It calls the optimization algorithm from pyswarm.
    """
    
    # Class-level lists for custom thicknesses
    int_thicklist = []
    long_thicklist = []
    
    # Class-level warning flags for optimization methods
    _flange_warning_logged = False
    _dimension_warning_logged = False
    _web_crippling_warning_logged = False

    def __init__(self):
        super(PlateGirderWelded, self).__init__()
        self.design_status = False
        # Instance-level warning flags
        self.flange_warning_logged = False
        self.dimension_warning_logged = False
        self.web_crippling_warning_logged = False
        
        # Attributes for optimization
        self._optimization_variable_list = []
        self._optimization_design_dictionary = {}
        self._optimization_is_thick_web = False
        self._optimization_is_symmetric = False


    ###############################################
    # Design Preference Functions Start
    ###############################################
    def tab_list(self):
        """
        Returns the list of tabs to be displayed in the Design Preferences window.
        """
        tabs = []

        t1 = (KEY_DISP_GIRDERSEC, TYPE_TAB_1, self.tab_girder_sec)
        tabs.append(t1)

        t5 = ("Optimisation", TYPE_TAB_2, self.optimization_tab_welded_plate_girder_design)
        tabs.append(t5)

        t1 = ("Stiffeners", TYPE_TAB_2, self.Stiffener_design)
        tabs.append(t1)

        t1 = ("Additional Girder Data", TYPE_TAB_2, self.girder_geometry)
        tabs.append(t1)

        t5 = ("Design", TYPE_TAB_2, self.design_values)
        tabs.append(t5)

        t6 = ("Deflection", TYPE_TAB_2, self.deflection_values)
        tabs.append(t6)

        return tabs

    def tab_value_changed(self):
        """
        Defines actions to be taken when values in Design Preferences tabs change.
        """
        change_tab = []

        t1 = (KEY_DISP_GIRDERSEC, [KEY_SEC_MATERIAL], [KEY_SEC_FU, KEY_SEC_FY], TYPE_TEXTBOX, self.get_fu_fy_I_section_plate_girder)
        change_tab.append(t1)

        t4 = (KEY_DISP_GIRDERSEC, ['Label_6', 'Label_7', 'Label_8', 'Label_9', 'Label_10', 'Label_11', KEY_SEC_FY],
              ['Label_12', 'Label_13', 'Label_14', 'Label_15', 'Label_16', 'Label_17', 'Label_18',
               'Label_19', 'Label_20', 'Label_21', 'Label_22', 'Label_23'], TYPE_TEXTBOX, self.Unsymm_I_Section_properties)
        change_tab.append(t4)

        t9 = ("Deflection", [KEY_STR_TYPE], [KEY_MEMBER_OPTIONS], TYPE_COMBOBOX, self.member_options_change)
        change_tab.append(t9)
        t9 = ("Deflection", [KEY_MEMBER_OPTIONS], [KEY_SUPPORTING_OPTIONS], TYPE_COMBOBOX, self.supp_options_change)
        change_tab.append(t9)
        t9 = ("Deflection", [KEY_STR_TYPE, KEY_DESIGN_LOAD, KEY_MEMBER_OPTIONS, KEY_SUPPORTING_OPTIONS], [KEY_MAX_DEFL], TYPE_TEXTBOX, self.max_defl_change)
        change_tab.append(t9)
        
        t10 = ("Stiffeners", [KEY_IntermediateStiffener_thickness], [KEY_IntermediateStiffener_thickness_val], TYPE_COMBOBOX, self.Int_stiffener_thickness_customized)
        change_tab.append(t10)
        t11 = ("Stiffeners", [KEY_LongitudnalStiffener_thickness], [KEY_LongitudnalStiffener_thickness_val], TYPE_COMBOBOX, self.Long_stiffener_thickness_customized)
        change_tab.append(t11)

        return change_tab

    def edit_tabs(self):
        """Not required for this module but empty list should be passed"""
        return []

    def input_dictionary_design_pref(self):
        """
        Selects values from Design Preferences to be saved to the design dictionary.
        """
        design_input = []

        t1 = (KEY_DISP_GIRDERSEC, TYPE_COMBOBOX, [KEY_SEC_MATERIAL])
        design_input.append(t1)
        t1 = (KEY_DISP_GIRDERSEC, TYPE_TEXTBOX, [KEY_SEC_FU, KEY_SEC_FY])
        design_input.append(t1)

        t2 = ("Optimisation", TYPE_TEXTBOX, [KEY_EFFECTIVE_AREA_PARA, KEY_LENGTH_OVERWRITE])
        design_input.append(t2)
        t2 = ("Optimisation", TYPE_COMBOBOX, [KEY_ALLOW_CLASS, KEY_LOAD])
        design_input.append(t2)

        t2 = ("Stiffeners", TYPE_COMBOBOX, [KEY_IntermediateStiffener, KEY_LongitudnalStiffener, KEY_IntermediateStiffener_thickness, KEY_LongitudnalStiffener_thickness])
        design_input.append(t2)
        t2 = ("Stiffeners", TYPE_TEXTBOX, [KEY_IntermediateStiffener_spacing])
        design_input.append(t2)
        t2 = ("Stiffeners", TYPE_COMBOBOX, [KEY_ShearBucklingOption, KEY_IntermediateStiffener_thickness_val, KEY_LongitudnalStiffener_thickness_val])
        design_input.append(t2)

        t2 = ("Additional Girder Data", TYPE_COMBOBOX, [KEY_IS_IT_SYMMETRIC])
        design_input.append(t2)

        t6 = ("Design", TYPE_COMBOBOX, [KEY_DP_DESIGN_METHOD])
        design_input.append(t6)

        t7 = ("Deflection", TYPE_COMBOBOX, [KEY_STR_TYPE, KEY_DESIGN_LOAD, KEY_MEMBER_OPTIONS, KEY_SUPPORTING_OPTIONS])
        design_input.append(t7)
        t7 = ("Deflection", TYPE_TEXTBOX, [KEY_MAX_DEFL])
        design_input.append(t7)

        return design_input

    def input_dictionary_without_design_pref(self):
        design_input = []
        t2 = (KEY_MATERIAL, [KEY_DP_DESIGN_METHOD], 'Input Dock')
        design_input.append(t2)
        t2 = (None, [KEY_ALLOW_CLASS, KEY_EFFECTIVE_AREA_PARA, KEY_LENGTH_OVERWRITE, KEY_LOAD, KEY_DP_DESIGN_METHOD, KEY_STR_TYPE, KEY_DESIGN_LOAD, KEY_MEMBER_OPTIONS, KEY_MAX_DEFL,
                     KEY_SUPPORTING_OPTIONS, KEY_ShearBucklingOption, KEY_IntermediateStiffener_spacing, KEY_IntermediateStiffener, KEY_LongitudnalStiffener, KEY_IntermediateStiffener_thickness_val, KEY_LongitudnalStiffener_thickness_val,
                     KEY_IntermediateStiffener_thickness, KEY_LongitudnalStiffener_thickness, KEY_IS_IT_SYMMETRIC], '')
        design_input.append(t2)
        return design_input

    def refresh_input_dock(self):
        add_buttons = []
        return add_buttons

    def get_values_for_design_pref(self, key, design_dictionary):
        """
        Provides default values for Design Preferences.
        """
        val = {
            KEY_ALLOW_CLASS: 'Yes',
            KEY_EFFECTIVE_AREA_PARA: '1.0',
            KEY_LENGTH_OVERWRITE: 'NA',
            KEY_LOAD: 'Normal',
            KEY_DP_DESIGN_METHOD: "Limit State Design",
            KEY_ShearBucklingOption: KEY_DISP_SB_Option[0],
            KEY_IS_IT_SYMMETRIC: 'Symmetrical',
            KEY_IntermediateStiffener_spacing: 'NA',
            KEY_IntermediateStiffener: 'No',
            KEY_IntermediateStiffener_thickness: 'All',
            KEY_LongitudnalStiffener: 'Yes and 1 stiffener',
            KEY_LongitudnalStiffener_thickness: 'All',
            KEY_STR_TYPE: 'Highway Bridge',
            KEY_DESIGN_LOAD: 'Live Load',
            KEY_MEMBER_OPTIONS: 'Simple Span',
            KEY_SUPPORTING_OPTIONS: 'NA',
            KEY_MAX_DEFL: 600,
            KEY_IntermediateStiffener_thickness_val: VALUES_STIFFENER_THICKNESS,
            KEY_LongitudnalStiffener_thickness_val: VALUES_STIFFENER_THICKNESS
        }[key]
        return val

    # --- Callbacks for Design Preferences ---

    def member_options_change(self):
        if self[0] == KEY_DISP_STR_TYP3:
            return {KEY_MEMBER_OPTIONS: VALUES_MEMBER_OPTIONS[1]}
        elif self[0] == KEY_DISP_STR_TYP4:
            return {KEY_MEMBER_OPTIONS: VALUES_MEMBER_OPTIONS[2]}
        else:
            return {KEY_MEMBER_OPTIONS: VALUES_MEMBER_OPTIONS[0]}

    def supp_options_change(self):
        if self[0] in ['Purlin and Girts', 'Simple span', 'Cantilever span']:
            return {KEY_SUPPORTING_OPTIONS: VALUES_SUPPORTING_OPTIONS_PSC}
        elif self[0] == 'Rafter Supporting':
            return {KEY_SUPPORTING_OPTIONS: VALUES_SUPPORTING_OPTIONS_RS}
        elif self[0] == 'Gantry':
            return {KEY_SUPPORTING_OPTIONS: VALUES_SUPPORTING_OPTIONS_GNT}
        elif self[0] in ['Floor and roof', 'Cantilever']:
            return {KEY_SUPPORTING_OPTIONS: VALUES_SUPPORTING_OPTIONS_FRC}
        else:
            return {KEY_SUPPORTING_OPTIONS: VALUES_SUPPORTING_OPTIONS_DEF}

    def max_defl_change(self):
        if self[0] in ['Highway Bridge', 'Railway Bridge']:
            if self[2] == 'Simple Span':
                if self[1] == 'Live load':
                    return {KEY_MAX_DEFL: VALUES_MAX_DEFL[0]}
                elif self[1] == 'Dead load':
                    return {KEY_MAX_DEFL: VALUES_MAX_DEFL[1]}
                else:
                    return {KEY_MAX_DEFL: 'NA'}
            else:
                if self[1] == 'Live load':
                    return {KEY_MAX_DEFL: VALUES_MAX_DEFL[2]}
                elif self[1] == 'Dead load':
                    return {KEY_MAX_DEFL: VALUES_MAX_DEFL[1]}
                else:
                    return {KEY_MAX_DEFL: 'NA'}
        elif self[0] == 'Other Building':
            if self[1] == 'Live load':
                if self[2] == 'Floor and roof':
                    if self[3] == 'Elements not susceptible to cracking':
                        return {KEY_MAX_DEFL: VALUES_MAX_DEFL[3]}
                    else:
                        return {KEY_MAX_DEFL: VALUES_MAX_DEFL[4]}
                else:
                    if self[3] == 'Elements not susceptible to cracking':
                        return {KEY_MAX_DEFL: VALUES_MAX_DEFL[5]}
                    else:
                        return {KEY_MAX_DEFL: VALUES_MAX_DEFL[6]}
            else:
                return {KEY_MAX_DEFL: 'NA'}
        else:
            if self[2] == 'Purlin and Girts' and self[1] == 'Live load':
                if self[3] == 'Elastic cladding':
                    return {KEY_MAX_DEFL: VALUES_MAX_DEFL[5]}
                else:
                    return {KEY_MAX_DEFL: VALUES_MAX_DEFL[6]}
            elif self[2] == 'Simple span' and self[1] == 'Live load':
                if self[3] == 'Elastic cladding':
                    return {KEY_MAX_DEFL: VALUES_MAX_DEFL[7]}
                else:
                    return {KEY_MAX_DEFL: VALUES_MAX_DEFL[3]}
            elif self[2] == 'Cantilever span' and self[1] == 'Live load':
                if self[3] == 'Elastic cladding':
                    return {KEY_MAX_DEFL: VALUES_MAX_DEFL[8]}
                else:
                    return {KEY_MAX_DEFL: VALUES_MAX_DEFL[5]}
            elif self[2] == 'Rafter Supporting' and self[1] == 'Live load':
                if self[3] == 'Profiled Metal sheeting':
                    return {KEY_MAX_DEFL: VALUES_MAX_DEFL[6]}
                else:
                    return {KEY_MAX_DEFL: VALUES_MAX_DEFL[7]}
            elif self[2] == 'Gantry' and self[1] == 'Live load':
                if self[1] == 'Crane Load(Manual operation)':
                    return {KEY_MAX_DEFL: VALUES_MAX_DEFL[9]}
                elif self[1] == 'Crane load(Electric operation up to 50t)':
                    return {KEY_MAX_DEFL: VALUES_MAX_DEFL[10]}
                else:
                    return {KEY_MAX_DEFL: VALUES_MAX_DEFL[11]}
            else:
                return {KEY_MAX_DEFL: 'NA'}

    def Int_stiffener_thickness_customized(self):
        """
        Launches the PopupDialog (from custom_widgets) for intermediate stiffeners.
        """
        selected_items = []
        if self[0] == 'All':
            return {KEY_IntermediateStiffener_thickness_val: VALUES_STIFFENER_THICKNESS}
        else:
            popup = PopupDialog()
            popup.listWidget.addItems(VALUES_STIFFENER_THICKNESS)  # Set available items
            if popup.exec_() == QDialog.Accepted:
                selected_items = popup.get_selected_items()
            PlateGirderWelded.int_thicklist = selected_items
            return {KEY_IntermediateStiffener_thickness_val: selected_items}

    def Long_stiffener_thickness_customized(self):
        """
        Launches the PopupDialog (from custom_widgets) for longitudinal stiffeners.
        """
        selected_items2 = []
        if self[0] == 'All':
            return {KEY_LongitudnalStiffener_thickness_val: VALUES_STIFFENER_THICKNESS}
        else:
            popup = PopupDialog()
            popup.listWidget.addItems(VALUES_STIFFENER_THICKNESS)  # Set available items
            if popup.exec_() == QDialog.Accepted:
                selected_items2 = popup.get_selected_items()
            PlateGirderWelded.long_thicklist = selected_items2
            return {KEY_LongitudnalStiffener_thickness_val: selected_items2}

    ####################################
    # Design Preference Functions End
    ####################################
    # Setting up logger and Input/Output Docks
    ####################################
    def module_name(self):
        return KEY_DISP_PLATE_GIRDER_WELDED

    @staticmethod
    def set_osdaglogger(key):
        """
        Set logger for Plate Girder Module.
        """
        global logger
        logger = logging.getLogger('Osdag')
        logger.setLevel(logging.DEBUG)
        
        # Avoid adding duplicate handlers
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter(fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            
            file_handler = logging.FileHandler('logging_text.log')
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

        if key is not None:
            # This logic might add duplicate OurLog handlers if called multiple times.
            # Consider adding a check to prevent duplicates if 'key' can change.
            handler = OurLog(key)
            formatter = logging.Formatter(fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
            handler.setFormatter(formatter)
            logger.addHandler(handler)

    def customized_input(self):
        """
        Defines custom input widgets for the Input Dock.
        """
        c_lst = []
        t1 = (KEY_TOP_FLANGE_THICKNESS_PG, self.plate_thick_customized)
        c_lst.append(t1)
        t2 = (KEY_BOTTOM_FLANGE_THICKNESS_PG, self.plate_thick_customized)
        c_lst.append(t2)
        t3 = (KEY_WEB_THICKNESS_PG, self.plate_thick_customized)
        c_lst.append(t3)
        return c_lst

    def input_values(self):
        """
        Defines the widgets and layout for the Input Dock.
        """
        self.module = KEY_DISP_PLATE_GIRDER_WELDED
        options_list = []

        t1 = (None, KEY_DISP_PG_SectionDetail, TYPE_TITLE, None, True, 'No Validator')
        options_list.append(t1)
        t1 = (KEY_MODULE, KEY_DISP_PLATE_GIRDER_WELDED, TYPE_MODULE, None, True, "No Validator")
        options_list.append(t1)
        t4 = (KEY_MATERIAL, KEY_DISP_MATERIAL, TYPE_COMBOBOX, VALUES_MATERIAL, True, 'No Validator')
        options_list.append(t4)
        t2 = (KEY_OVERALL_DEPTH_PG_TYPE, KEY_DISP_OVERALL_DEPTH_PG_TYPE, TYPE_COMBOBOX, VALUES_DEPTH_PG, True, 'No Validator')
        options_list.append(t2)
        t33 = (KEY_OVERALL_DEPTH_PG, KEY_DISP_OVERALL_DEPTH_PG, TYPE_TEXTBOX, None, True, 'Int Validator')
        options_list.append(t33)
        t4 = (KEY_WEB_THICKNESS_PG, KEY_DISP_WEB_THICKNESS_PG, TYPE_COMBOBOX_CUSTOMIZED, VALUES_PLATETHK, True, 'Int Validator')
        options_list.append(t4)
        t2 = (KEY_TOP_Bflange_PG, KEY_DISP_TOP_Bflange_PG, TYPE_TEXTBOX, None, True, 'Int Validator')
        options_list.append(t2)
        t4 = (KEY_TOP_FLANGE_THICKNESS_PG, KEY_DISP_TOP_FLANGE_THICKNESS_PG, TYPE_COMBOBOX_CUSTOMIZED, VALUES_PLATETHK, True, 'Int Validator')
        options_list.append(t4)
        t22 = (KEY_BOTTOM_Bflange_PG, KEY_DISP_BOTTOM_Bflange_PG, TYPE_TEXTBOX, None, True, 'Int Validator')
        options_list.append(t22)
        t4 = (KEY_BOTTOM_FLANGE_THICKNESS_PG, KEY_DISP_BOTTOM_FLANGE_THICKNESS_PG, TYPE_COMBOBOX_CUSTOMIZED, VALUES_PLATETHK, True, 'No Validator')
        options_list.append(t4)
        t2 = (KEY_LENGTH, KEY_DISP_LENGTH, TYPE_TEXTBOX, None, True, 'No Validator')
        options_list.append(t2)

        t1 = (None, KEY_DISP_SECTION_DATA_PG, TYPE_TITLE, None, True, 'No Validator')
        options_list.append(t1)
        t2 = (KEY_DESIGN_TYPE_FLEXURE, KEY_BEAM_SUPP_TYPE, TYPE_COMBOBOX, VALUES_SUPP_TYPE_temp, True, "No Validator")
        options_list.append(t2)
        t5 = (KEY_SUPPORT_WIDTH, KEY_DISP_SUPPORT_WIDTH, TYPE_TEXTBOX, None, True, 'Int Validator')
        options_list.append(t5)
        t4 = (KEY_WEB_PHILOSOPHY, KEY_DISP_WEB_PHILOSOPHY, TYPE_COMBOBOX, WEB_PHILOSOPHY_list, True, 'No Validator')
        options_list.append(t4)
        t10 = (KEY_TORSIONAL_RES, DISP_TORSIONAL_RES, TYPE_COMBOBOX, Torsion_Restraint_list, True, 'No Validator')
        options_list.append(t10)
        t11 = (KEY_WARPING_RES, DISP_WARPING_RES, TYPE_COMBOBOX, Warping_Restraint_list, True, 'No Validator')
        options_list.append(t11)

        t7 = (None, KEY_LOADING, TYPE_TITLE, None, True, 'No Validator')
        options_list.append(t7)
        t8 = (KEY_MOMENT, KEY_DISP_MOMENT, TYPE_TEXTBOX, None, True, 'No Validator')
        options_list.append(t8)
        t8 = (KEY_SHEAR, KEY_DISP_SHEAR, TYPE_TEXTBOX, None, True, 'No Validator')
        options_list.append(t8)
        t8 = (KEY_BENDING_MOMENT_SHAPE, KEY_DISP_BENDING_MOMENT_SHAPE, TYPE_COMBOBOX, Bending_moment_shape_list, True, 'No Validator')
        options_list.append(t8)

        return options_list

    # --- Callbacks for Input Dock ---

    def fn_torsion_warping(self):
        if self[0] == Torsion_Restraint1:
            return Warping_Restraint_list
        elif self[0] == Torsion_Restraint2:
            return [Warping_Restraint5]
        else:
            return [Warping_Restraint5]

    def axis_bending_change(self):
        if self[0] == KEY_DISP_DESIGN_TYPE_FLEXURE:
            return ['NA']
        else:
            return VALUES_BENDING_TYPE

    def fn_conn_image(self):
        img = self[0]
        if img == Bending_moment_shape_list[0]:
            return VALUES_IMAGE_PLATEGIRDER[0]
        elif img == Bending_moment_shape_list[1]:
            return VALUES_IMAGE_PLATEGIRDER[1]
        elif img == Bending_moment_shape_list[2]:
            return VALUES_IMAGE_PLATEGIRDER[2]
        elif img == Bending_moment_shape_list[3]:
            return VALUES_IMAGE_PLATEGIRDER[3]
        else:
            return VALUES_IMAGE_PLATEGIRDER[4]

    def customized_dimensions(self):
        return KEY_DISP_OVERALL_DEPTH_PG if self[0] == "Customized" else ''

    def customized_dimensions_1(self):
        return KEY_DISP_TOP_Bflange_PG if self[0] == "Customized" else ''

    def customized_dimensions_2(self):
        return KEY_DISP_BOTTOM_Bflange_PG if self[0] == "Customized" else ''

    def customized_dims(self):
        return True if self[0] == "Customized" else False

    def customized_options(self):
        return VALUES_PLATETHK if self[0] == "Customized" else VALUES_OPT

    def customized_dimensions_cst(self):
        return KEY_OVERALL_DEPTH_PG_CST if self[0] == "Optimized" else ''

    def customized_dims_cst(self):
        return True if self[0] == "Optimized" else False

    def pop_up_bounds(self):
        """
        Launches the RangeInputDialog (from custom_widgets).
        """
        if self[0] == "Bound Values":
            dialog = RangeInputDialog()
            if dialog.exec_() == QDialog.Accepted:
                return str(dialog.get_values())

    def input_value_changed(self):
        """
        Defines actions to be taken when values in the Input Dock change.
        """
        lst = []
        t3 = ([KEY_TORSIONAL_RES], KEY_WARPING_RES, TYPE_COMBOBOX, self.fn_torsion_warping)
        lst.append(t3)
        t44 = ([KEY_OVERALL_DEPTH_PG_TYPE], KEY_OVERALL_DEPTH_PG, TYPE_LABEL, self.customized_dimensions)
        lst.append(t44)
        t45 = ([KEY_OVERALL_DEPTH_PG_TYPE], KEY_OVERALL_DEPTH_PG, TYPE_TEXTBOX, self.customized_dims)
        lst.append(t45)
        t2 = ([KEY_OVERALL_DEPTH_PG_TYPE], KEY_TOP_Bflange_PG, TYPE_LABEL, self.customized_dimensions_1)
        lst.append(t2)
        t3 = ([KEY_OVERALL_DEPTH_PG_TYPE], KEY_TOP_Bflange_PG, TYPE_TEXTBOX, self.customized_dims)
        lst.append(t3)
        t23 = ([KEY_OVERALL_DEPTH_PG_TYPE], KEY_BOTTOM_Bflange_PG, TYPE_LABEL, self.customized_dimensions_2)
        lst.append(t23)
        t24 = ([KEY_OVERALL_DEPTH_PG_TYPE], KEY_BOTTOM_Bflange_PG, TYPE_TEXTBOX, self.customized_dims)
        lst.append(t24)
        t3 = ([KEY_MATERIAL], KEY_MATERIAL, TYPE_CUSTOM_MATERIAL, self.new_material)
        lst.append(t3)
        
        # Output dock visibility modifiers
        for key in [KEY_T_constatnt, KEY_W_constatnt, KEY_Elastic_CM]:
            lst.append(([KEY_DESIGN_TYPE_FLEXURE], key, TYPE_OUT_LABEL, self.output_modifier))
            lst.append(([KEY_DESIGN_TYPE_FLEXURE], key, TYPE_OUT_DOCK, self.output_modifier))
            
        for key in [KEY_IntermediateStiffener_thickness, KEY_LongitudnalStiffener_thickness,
                    KEY_IntermediateStiffener_spacing, KEY_LongitudnalStiffener_numbers,
                    KEY_LongitudinalStiffener1_pos, KEY_LongitudinalStiffener2_pos]:
            lst.append(([KEY_WEB_PHILOSOPHY], key, TYPE_OUT_LABEL, self.output_modifier2))
            lst.append(([KEY_WEB_PHILOSOPHY], key, TYPE_OUT_DOCK, self.output_modifier2))

        return lst

    def warning_majorbending(self):
        return True if self[0] == VALUES_SUPP_TYPE_temp[2] else False

    def output_modifier(self):
        return False if self[0] == VALUES_SUPP_TYPE_temp[2] else True
        
    def output_modifier_long_stiffener(self):
        return False if self[0] == 'Thin we' else True
        
    def output_modifier2(self):
        # This seems to be the wrong logic, but keeping it as it was.
        # It hides stiffener info if 'Thin Web with ITS' is selected.
        return False if self[0] == 'Thin Web with ITS' else True

    def output_values(self, flag):
        """
        Defines the widgets and layout for the Output Dock.
        """
        out_list = []
        t0 = (None, DISP_TITLE_STRUT_SECTION, TYPE_TITLE, None, True)
        out_list.append(t0)

        t1 = (KEY_TITLE_OPTIMUM_DESIGNATION, KEY_DISP_TITLE_OPTIMUM_DESIGNATION, TYPE_TEXTBOX,
              self.result_designation if flag else '', True)
        out_list.append(t1)
        t2 = (KEY_OPTIMUM_UR_COMPRESSION, KEY_DISP_OPTIMUM_UR_COMPRESSION, TYPE_TEXTBOX, round(self.result_UR, 3) if flag else '', True)
        out_list.append(t2)
        t3 = (KEY_OPTIMUM_SC, KEY_DISP_OPTIMUM_SC, TYPE_TEXTBOX, self.section_classification_val if flag else '', True)
        out_list.append(t3)
        t4 = (KEY_betab_constatnt, KEY_DISP_betab_constatnt, TYPE_TEXTBOX,
              self.betab if flag else '', True)
        out_list.append(t4)
        t5 = (KEY_EFF_SEC_AREA, KEY_DISP_EFF_SEC_AREA, TYPE_TEXTBOX, self.effectivearea if flag else '',
              True)
        out_list.append(t5)

        t10 = (KEY_IntermediateStiffener_thickness, KEY_DISP_IntermediateStiffener_thickness, TYPE_TEXTBOX,
               self.intstiffener_thk if flag else '', True)
        out_list.append(t10)
        t10 = (KEY_IntermediateStiffener_spacing, KEY_DISP_IntermediateStiffener_spacing, TYPE_TEXTBOX,
               self.intstiffener_spacing if flag else '', True)
        out_list.append(t10)
        t1 = (KEY_LongitudnalStiffener_thickness, KEY_DISP_LongitudnalStiffener_thickness, TYPE_TEXTBOX,
              self.longstiffener_thk if flag else '', True)
        out_list.append(t1)
        t1 = (KEY_LongitudnalStiffener_numbers, KEY_DISP_LongitudnalStiffener_numbers, TYPE_TEXTBOX, self.longstiffener_no if flag else '', True)
        out_list.append(t1)
        t2 = (KEY_EndpanelStiffener_thickness, KEY_DISP_EndpanelStiffener_thickness, TYPE_TEXTBOX, self.end_panel_stiffener_thickness if flag else '', True)
        out_list.append(t2)
        t1 = (KEY_MOMENT_STRENGTH, KEY_DISP_MOMENT, TYPE_TEXTBOX,
              self.design_moment if flag else '', True)
        out_list.append(t1)

        t1 = (KEY_WeldWebtoflange, KEY_DISP_WeldWebtoflange, TYPE_TEXTBOX,
              max(self.atop, self.abot) if flag else '', True)
        out_list.append(t1)
        t1 = (KEY_WeldStiffenertoweb, KEY_DISP_WeldStiffenertoweb, TYPE_TEXTBOX,
              self.weld_stiff if flag else '', True)
        out_list.append(t1)

        t2 = (KEY_T_constatnt, KEY_DISP_T_constatnt, TYPE_TEXTBOX,
              self.torsion_cnst if flag else '', False)
        out_list.append(t2)
        t2 = (KEY_W_constatnt, KEY_DISP_W_constatnt, TYPE_TEXTBOX, self.warping_cnst if flag else '', False)
        out_list.append(t2)
        t2 = (KEY_LongitudinalStiffener1_pos, KEY_DISP_LongitudinalStiffener1_pos, TYPE_TEXTBOX, self.x1 if flag else '', True)
        out_list.append(t2)
        t2 = (KEY_LongitudinalStiffener2_pos, KEY_DISP_LongitudinalStiffener2_pos, TYPE_TEXTBOX, self.x2 if flag else '', True)
        out_list.append(t2)
        t2 = (KEY_Elastic_CM, KEY_DISP_Elastic_CM, TYPE_TEXTBOX, self.critical_moment if flag else '', False)
        out_list.append(t2)

        return out_list

    def spacing(self, status):
        """
        Defines additional output values (legacy or for specific tabs).
        """
        spacing = []
        t2 = (KEY_T_constatnt, KEY_DISP_T_constatnt, TYPE_TEXTBOX,
              self.result_tc if status else '', False)
        spacing.append(t2)
        t2 = (KEY_W_constatnt, KEY_DISP_W_constatnt, TYPE_TEXTBOX, self.result_wc if status else '', False)
        spacing.append(t2)
        t2 = (KEY_IMPERFECTION_FACTOR_LTB, KEY_DISP_IMPERFECTION_FACTOR, TYPE_TEXTBOX, self.result_IF_lt if status else '',
              False)
        spacing.append(t2)
        t2 = (KEY_SR_FACTOR_LTB, KEY_DISP_SR_FACTOR, TYPE_TEXTBOX, self.result_srf_lt if status else '', False)
        spacing.append(t2)
        t2 = (KEY_NON_DIM_ESR_LTB, KEY_DISP_NON_DIM_ESR, TYPE_TEXTBOX, self.result_nd_esr_lt if status else '', False)
        spacing.append(t2)
        t1 = (KEY_DESIGN_STRENGTH_COMPRESSION, KEY_DISP_COMP_STRESS, TYPE_TEXTBOX,
              self.result_fcd__lt if status else '', False)
        spacing.append(t1)
        t2 = (KEY_Elastic_CM, KEY_DISP_Elastic_CM, TYPE_TEXTBOX, self.result_mcr if status else '', False)
        spacing.append(t2)
        return spacing

    ####################################
    # Main Design Orchestration
    ####################################

    def func_for_validation(self, main_arg, design_dictionary):
        """
        Validates all inputs before starting the design.
        
        NOTE: Accepts a redundant 'main_arg' to match a bug in the Osdag
        framework's calling convention in ui_template.py.
        """
        all_errors = []
        self.design_status = False
        flag = False
        self.output_values(flag) # This call is now correct
        
        flags = {KEY_LENGTH: False, KEY_SHEAR: False, KEY_MOMENT: False}
        option_list = self.input_values()
        missing_fields_list = []

        for option in option_list:
            key, display_name, ui_type, _, required, _ = option
            if not required:
                continue

            value = design_dictionary.get(key)

            if ui_type == TYPE_TEXTBOX:
                if not value:
                    # Special check for Optimized mode
                    if design_dictionary.get(KEY_OVERALL_DEPTH_PG_TYPE) == 'Optimized' and \
                       key in [KEY_OVERALL_DEPTH_PG, KEY_TOP_Bflange_PG, KEY_BOTTOM_Bflange_PG]:
                        pass # These fields are allowed to be empty in Optimized mode
                    else:
                        missing_fields_list.append(display_name)
                    continue
                
                try:
                    float_val = float(value)
                    if float_val <= 0.0:
                        all_errors.append(f"Input value for '{display_name}' must be greater than zero.")
                    if key in flags:
                        flags[key] = True
                except ValueError:
                    all_errors.append(f"Input value for '{display_name}' is not a valid number.")
            
            elif ui_type in [TYPE_COMBOBOX, TYPE_COMBOBOX_CUSTOMIZED]:
                 if not value or value == "Select" or value == []:
                     missing_fields_list.append(display_name)

        if len(missing_fields_list) > 0:
            error = self.generate_missing_fields_error_string(missing_fields_list)
            all_errors.append(error)
        
        if not all(flags.values()):
             all_errors.append("Length, Shear, and Moment must be provided and be > 0.")

        if not all_errors:
            # All checks passed, call set_input_values which triggers the design
            self.set_input_values(design_dictionary)
        else:
            return all_errors

    def set_input_values(self, design_dictionary):
        """
        Sets up 'self' attributes from the validated design_dictionary
        and then calls the main design/optimization method.
        """
        self.module = design_dictionary[KEY_MODULE]
        self.mainmodule = 'PLATE GIRDER'
        self.design_type = design_dictionary[KEY_OVERALL_DEPTH_PG_TYPE]
        self.section_class = None

        if self.design_type == 'Optimized':
            self.total_depth = 1 # Placeholder
            self.web_thickness_list = design_dictionary[KEY_WEB_THICKNESS_PG]
            self.top_flange_width = 1 # Placeholder
            self.top_flange_thickness_list = design_dictionary[KEY_TOP_FLANGE_THICKNESS_PG]
            self.bottom_flange_width = 1 # Placeholder
            self.bottom_flange_thickness_list = design_dictionary[KEY_BOTTOM_FLANGE_THICKNESS_PG]
            
            # Initialize with first value for material property check
            self.web_thickness = float(design_dictionary[KEY_WEB_THICKNESS_PG][0])
            self.top_flange_thickness = float(design_dictionary[KEY_TOP_FLANGE_THICKNESS_PG][0])
            self.bottom_flange_thickness = float(design_dictionary[KEY_BOTTOM_FLANGE_THICKNESS_PG][0])
        else:
            self.total_depth = float(design_dictionary[KEY_OVERALL_DEPTH_PG])
            self.web_thickness_list = design_dictionary[KEY_WEB_THICKNESS_PG]
            self.web_thickness = float(design_dictionary[KEY_WEB_THICKNESS_PG][0])
            self.top_flange_width = float(design_dictionary[KEY_TOP_Bflange_PG])
            self.top_flange_thickness_list = design_dictionary[KEY_TOP_FLANGE_THICKNESS_PG]
            self.top_flange_thickness = float(design_dictionary[KEY_TOP_FLANGE_THICKNESS_PG][0])
            self.bottom_flange_width = float(design_dictionary[KEY_BOTTOM_Bflange_PG])
            self.bottom_flange_thickness_list = design_dictionary[KEY_BOTTOM_FLANGE_THICKNESS_PG]
            self.bottom_flange_thickness = float(design_dictionary[KEY_BOTTOM_FLANGE_THICKNESS_PG][0])

        thickness_for_mat = max(self.web_thickness, self.top_flange_thickness, self.bottom_flange_thickness)
        self.material = Material(design_dictionary[KEY_MATERIAL], thickness_for_mat)
        
        if self.total_depth > self.top_flange_thickness + self.bottom_flange_thickness:
            self.eff_depth = self.total_depth - self.top_flange_thickness - self.bottom_flange_thickness
        else:
            self.eff_depth = 0 # Invalid
            
        self.IntStiffnerwidth = min(self.top_flange_width, self.bottom_flange_width) - self.web_thickness / 2 - 10
        self.eff_width_longitudnal = min(self.top_flange_width, self.bottom_flange_width) - self.web_thickness / 2 - 10
        
        # Set custom thickness lists
        if design_dictionary[KEY_IntermediateStiffener_thickness] == 'Customized':
            design_dictionary[KEY_IntermediateStiffener_thickness_val] = PlateGirderWelded.int_thicklist
        else:
            design_dictionary[KEY_IntermediateStiffener_thickness_val] = VALUES_STIFFENER_THICKNESS
        self.int_thickness_list = design_dictionary[KEY_IntermediateStiffener_thickness_val]

        if design_dictionary[KEY_LongitudnalStiffener_thickness] == 'Customized':
            design_dictionary[KEY_LongitudnalStiffener_thickness_val] = PlateGirderWelded.long_thicklist
        else:
            design_dictionary[KEY_LongitudnalStiffener_thickness_val] = VALUES_STIFFENER_THICKNESS
        self.long_thickness_list = design_dictionary[KEY_LongitudnalStiffener_thickness_val]
        
        self.deflection_criteria = design_dictionary[KEY_MAX_DEFL]
        self.support_condition = 'Simply Supported'
        self.loading_case = design_dictionary[KEY_BENDING_MOMENT_SHAPE]
        self.shear_type = None
        self.support_type = design_dictionary[KEY_DESIGN_TYPE_FLEXURE]
        self.loading_condition = design_dictionary[KEY_LOAD]
        self.torsional_res = design_dictionary[KEY_TORSIONAL_RES]
        self.warping = design_dictionary[KEY_WARPING_RES]
        self.length = float(design_dictionary[KEY_LENGTH])
        self.effective_length = None
        self.allow_class = design_dictionary[KEY_ALLOW_CLASS]
        self.beta_b_lt = None
        self.web_philosophy = design_dictionary[KEY_WEB_PHILOSOPHY]
        self.epsilon = math.sqrt(250 / (self.material.fy)) # fy is in MPa
        self.b1 = float(design_dictionary[KEY_SUPPORT_WIDTH])
        self.c = design_dictionary[KEY_IntermediateStiffener_spacing]
        self.Is = None
        self.IntStiffThickness = float(self.int_thickness_list[0]) if self.int_thickness_list else 6.0
        self.LongStiffThickness = float(self.long_thickness_list[0]) if self.long_thickness_list else 6.0
        self.x1 = 0
        self.x2 = 0
        self.V_cr = None
        self.V_d = None
        self.V_tf = None
        self.long_Stiffner = design_dictionary[KEY_LongitudnalStiffener]
        self.load = Load(shear_force=design_dictionary[KEY_SHEAR], axial_force="", moment=design_dictionary[KEY_MOMENT], unit_kNm=True)
        self.alpha_lt = 0.49 # for welded sections
        self.phi_lt = None
        self.gamma_m0 = IS800_2007.cl_5_4_1_Table_5["gamma_m0"]["yielding"]
        self.X_lt = None
        self.fbd_lt = None
        self.Md = None
        self.lefactor = 0.7 # Per 8.7.2.4 for intermediate stiffeners
        self.M_cr = None
        self.F_q = None
        self.Critical_buckling_load = None
        self.shear_ratio = 0
        self.endshear_ratio = 0
        self.moment_ratio = 0
        self.deflection_ratio = 0
        self.web_buckling_ratio = 0
        self.It = None
        self.Iw = None
        self.torsion_cnst = None
        self.warping_cnst = None
        self.critical_moment = None
        self.fcd = None
        self.end_stiffthickness = 0
        self.stiffener_type = None
        self.end_panel_stiffener_thickness = None
        self.end_stiffwidth = min(self.top_flange_width, self.bottom_flange_width) / 2 - self.web_thickness / 2 - 10
        if self.end_stiffwidth <= 0: self.end_stiffwidth = 50 # Fallback
        
        self.design_status = False

        # --- Trigger Design ---
        if self.design_type == 'Optimized':
            is_thick_web = self.web_philosophy == 'Thick Web without ITS'
            is_symmetric = design_dictionary[KEY_IS_IT_SYMMETRIC] == 'Symmetrical'
            
            # Call the inherited optimization method
            self.optimized_method(design_dictionary, is_thick_web, is_symmetric)
        else:
            # Call the inherited design check method
            self.design_check(design_dictionary, silent=False)

    def design_check(self, design_dictionary, silent=False):
        """
        Orchestrates the design checks for a "Customized" girder.
        All calculation methods are inherited from PlateGirderLogic.
        """
        self.design_flag = False
        self.design_flag2 = False
        self.shearflag1 = False
        self.shearflag2 = False
        self.shearflag3 = False
        self.shearchecks = False
        self.momentchecks = False
        self.defl_check = False
        self.long_check = False
        
        # 0. Check for valid geometry
        if self.eff_depth <= 0:
            if not silent: logger.error("Invalid Geometry: Total Depth is less than combined flange thickness.")
            return False
        
        # 1. Section Classification
        self.design_flag = self.section_classification(design_dictionary)
        if not self.design_flag:
            if not silent: logger.error("Slender section not allowed for 'Thick Web' philosophy.")
            return False
            
        # Log flange and dimension warnings (only if not silent)
        if not silent:
            if self.top_flange_thickness > 0:
                min_b_tf = 7.4 * self.epsilon
                if (self.top_flange_width / self.top_flange_thickness) < min_b_tf:
                    logger.warning("Top flange b/tf ratio is low, flanges may be too thick.")
                if (self.bottom_flange_width / self.bottom_flange_thickness) < min_b_tf:
                    logger.warning("Bottom flange b/tf ratio is low, flanges may be too thick.")
            if self.bottom_flange_width < self.top_flange_width:
                logger.warning("Bottom flange width is less than top flange width.")
            if self.bottom_flange_thickness < self.top_flange_thickness:
                logger.warning("Bottom flange thickness is less than top flange thickness.")

        # 2. Beta Value
        self.beta_value(design_dictionary, self.section_class)

        # 3. Strength Checks (Shear & Moment)
        if self.web_philosophy == 'Thick Web without ITS':
            self.design_flag2 = self.min_web_thickness_thick_web(self, self.eff_depth, self.web_thickness, self.epsilon, "no_stiffener", 0, silent=silent)
            if self.design_flag2:
                self.shearflag1 = self.shear_capacity_laterally_supported_thick_web(self, self.material.fy, self.gamma_m0, self.total_depth, self.web_thickness, self.top_flange_thickness, self.bottom_flange_thickness)
                self.shearflag2 = self.web_buckling_laterally_supported_thick_web(self, self.material.fy, self.gamma_m0, self.total_depth, self.web_thickness, self.top_flange_thickness, self.bottom_flange_thickness, self.material.modulus_of_elasticity, self.b1)
                self.shearflag3 = self.web_crippling_laterally_supported_thick_web(self, self.material.fy, self.gamma_m0, self.web_thickness, self.top_flange_thickness, self.b1, silent=silent)
                self.shearchecks = self.shearflag1 and self.shearflag2 and self.shearflag3
                
                if self.support_type == 'Major Laterally Supported':
                    self.momentchecks = self.moment_capacity_laterally_supported(self, self.load.shear_force, self.plast_sec_mod_z, self.elast_sec_mod_z, self.material.fy, self.gamma_m0, self.total_depth, self.web_thickness, self.top_flange_thickness, self.bottom_flange_thickness, self.section_class)
                else:
                    self.momentchecks = self.moment_capacity_laterally_unsupported(self, self.material.modulus_of_elasticity, self.effective_length, self.total_depth, self.top_flange_thickness, self.bottom_flange_thickness, self.top_flange_width, self.bottom_flange_width, self.web_thickness, self.loading_case, self.gamma_m0, self.material.fy, self.load.shear_force, silent=silent)
            else:
                if not silent: logger.error("Web thickness is insufficient for 'Thick Web' philosophy.")
                return False

        else: # Thin Web
            self.shear_ratio = 0
            if self.long_Stiffner == 'Yes and 1 stiffener': self.stiffener_type = "transverse_and_one_longitudinal_compression"
            elif self.long_Stiffner == 'Yes and 2 stiffeners': self.stiffener_type = "transverse_and_two_longitudinal_neutral"
            else: self.stiffener_type = "transverse_only"
            
            if self.c == 'NA':
                if not silent: logger.error("Stiffener spacing (c) not provided for 'Thin Web' design.")
                return False
            
            self.c = float(self.c)
            if self.c <= 0:
                if not silent: logger.error("Stiffener spacing (c) must be greater than 0.")
                return False
                
            if self.stiffener_type != "transverse_only":
                self.long_check = self.design_longitudinal_stiffeners(self, self.eff_depth, self.web_thickness, self.c, self.epsilon, (self.stiffener_type == "transverse_and_two_longitudinal_neutral"), silent=silent)
                if not self.long_check and not silent: logger.error("Longitudinal Stiffener Check failed")
            
            self.design_flag2 = self.min_web_thickness_thick_web(self, self.eff_depth, self.web_thickness, self.epsilon, self.stiffener_type, self.c, silent=silent)
            if self.design_flag2:
                if design_dictionary[KEY_ShearBucklingOption] == 'Simple Post Critical':
                    self.shearflag1 = self.shear_buckling_check_simple_postcritical(self, self.eff_depth, self.total_depth, self.top_flange_thickness, self.bottom_flange_thickness, self.web_thickness, self.load.shear_force, self.c)
                    self.shearflag2 = self.shear_buckling_check_intermediate_stiffener(self, self.eff_depth, self.web_thickness, self.c, self.epsilon, self.IntStiffThickness, self.IntStiffnerwidth, self.load.shear_force, self.gamma_m0, self.material.fy, self.material.modulus_of_elasticity)
                else: # Tension Field
                    self.shearflag1 = self.shear_buckling_check_tension_field(self, self.eff_depth, self.total_depth, self.top_flange_thickness, self.bottom_flange_thickness, self.web_thickness, self.c)
                    self.shearflag2 = self.tension_field_intermediate_stiffener(self, self.eff_depth, self.web_thickness, self.c, self.epsilon, self.IntStiffThickness, self.IntStiffnerwidth, self.load.shear_force, self.gamma_m0, self.material.fy, self.material.modulus_of_elasticity)

                self.shearchecks = self.shearflag1 and self.shearflag2
                
                if self.support_type == 'Major Laterally Supported':
                    self.momentchecks = self.moment_capacity_laterally_supported(self, self.load.shear_force, self.plast_sec_mod_z, self.elast_sec_mod_z, self.material.fy, self.gamma_m0, self.total_depth, self.web_thickness, self.top_flange_thickness, self.bottom_flange_thickness, self.section_class)
                else:
                    self.momentchecks = self.moment_capacity_laterally_unsupported(self, self.material.modulus_of_elasticity, self.effective_length, self.total_depth, self.top_flange_thickness, self.bottom_flange_thickness, self.top_flange_width, self.bottom_flange_width, self.web_thickness, self.loading_case, self.gamma_m0, self.material.fy, self.load.shear_force, silent=silent)
            else:
                if not silent: logger.error("Web thickness is insufficient for 'Thin Web' philosophy (with stiffeners).")
                return False

        # 4. End Panel Stiffener
        if not self.end_panel_stiffener_calc(self, self.top_flange_width, self.bottom_flange_width, self.web_thickness, self.material.fy, self.gamma_m0, self.eff_depth, self.top_flange_thickness, self.total_depth, self.effective_length, self.bottom_flange_thickness, self.material.modulus_of_elasticity, self.epsilon, self.c, silent=silent):
            if not silent: logger.error("End Panel Stiffener Check failed")
            self.shearchecks = False # End panel failure is a shear failure
        elif not silent:
             logger.info("End Panel Stiffener Check passed")
             
        # 5. Deflection check
        self.defl_check = self.evaluate_deflection_kNm_mm(self, self.load.moment, self.length, self.material.modulus_of_elasticity, self.loading_case, self.deflection_criteria, silent=silent)
        if not self.defl_check and not silent: logger.error("Deflection Check failed")

        # 6. Final Result
        all_checks_passed = self.momentchecks and self.shearchecks and self.defl_check
        
        if not silent:
            if all_checks_passed:
                logger.info("Design is SAFE.")
                self.final_format(design_dictionary)
                self.design_status = True
            else:
                logger.error("Design FAILED. Check logs for details.")
                self.final_format(design_dictionary) # Format with failing values
                self.design_status = False
        
        return all_checks_passed


    def _objective_function_with_penalty(self, particle):
        """
        Internal helper for PSO.
        Calculates mass and adds a penalty for failed checks.
        """
        # 1. Assign particle dimensions to self
        sec = self.assign_particle_to_section(particle, self._optimization_variable_list)
        
        # 2. Calculate mass (the objective)
        area = (sec.bf_top * sec.tf_top +
                sec.bf_bot * sec.tf_bot +
                sec.tw * (sec.D - sec.tf_top - sec.tf_bot))
        mass = area * self.length * 7.85e-6  # kg

        if not self._optimization_is_thick_web:
            if sec.c > 0: # Avoid division by zero
                n_stiff = max(self.length / sec.c - 1, 0)
                width_stiff = min(sec.bf_top, sec.bf_bot) - sec.tw / 2 - 10
                height_stiff = sec.D - sec.tf_top - sec.tf_bot
                if width_stiff > 0 and height_stiff > 0:
                    vol_stiff = n_stiff * 2 * width_stiff * sec.t_stiff * height_stiff
                    mass_stiff = vol_stiff * 7.85e-6
                    mass += mass_stiff
            
        # 3. Run all checks silently
        # We must pass a copy of the dictionary, as design_check can modify it
        design_dict_copy = self._optimization_design_dictionary.copy()
        is_safe = self.design_check(design_dict_copy, silent=True)

        # 4. Apply Penalty
        penalty = 0.0
        if not is_safe:
            # Add a penalty based on the worst-failing ratio
            # This guides the optimizer back to a safe region
            max_ratio = max(self.moment_ratio, self.shear_ratio, self.deflection_ratio)
            if max_ratio > 1.0:
                penalty = (max_ratio - 1.0) * 1e6 # High penalty for failing
            else:
                penalty = 1e6 # High fixed penalty if a check failed but ratios are < 1 (e.g., thickness)

        return mass + penalty


    def optimized_method(self, design_dictionary, is_thick_web, is_symmetric):
        """
        Orchestrates the Particle Swarm Optimization using pyswarm.
        """
        logger.info("Starting PSO optimization...")
        
        # Reset warning flags
        self.flange_warning_logged = False
        self.dimension_warning_logged = False
        self.web_crippling_warning_logged = False
        PlateGirderWelded._web_crippling_warning_logged = False

        # Store context for the objective function
        self._optimization_variable_list = self.build_variable_structure(is_thick_web, is_symmetric)
        self._optimization_design_dictionary = design_dictionary
        self._optimization_is_thick_web = is_thick_web
        self._optimization_is_symmetric = is_symmetric
        
        lb, ub = self.get_bounds(self._optimization_variable_list)
        
        # Generate a feasible first particle
        initial_particle = self.generate_first_particle(
            self, self.length, self.load.moment, self.material.fy,
            is_thick_web, is_symmetric
        )
        # Ensure the initial particle is within the bounds
        initial_particle = np.clip(initial_particle, lb, ub)
        
        # Use a list to pass the initial swarm position
        initial_swarm = [initial_particle]
        
        # Run the PSO from pyswarm
        # We use a simple penalty function, so no `f_ieqcons` is needed.
        best_pos, best_cost = pso(
            self._objective_function_with_penalty,
            lb,
            ub,
            swarmsize=50,
            maxiter=50,
            minstep=1e-6,
            minfunc=1e-6,
            debug=False,
            particle_output=True, # pyswarm needs this to accept initial swarm
            pso_locations = initial_swarm # Pass our guessed particle
        )
        
        logger.info(f"PSO calculation successfully completed with best cost: {best_cost}")
        
        # Apply optimized values and round up to standard plates/dimensions
        best_design_var = dict(zip(self._optimization_variable_list, best_pos))

        def ceil_to_nearest(x, multiple):
            return float(math.ceil(x / multiple) * multiple)

        if is_symmetric:
            opt_tf = float(best_design_var['tf'])
            self.bottom_flange_thickness = self.top_flange_thickness = next((t for t in self.top_flange_thickness_list if float(t) >= opt_tf), self.top_flange_thickness_list[-1])
            
            opt_tw = float(best_design_var['tw'])
            self.web_thickness = next((t for t in self.web_thickness_list if float(t) >= opt_tw), self.web_thickness_list[-1])

            self.top_flange_width = self.bottom_flange_width = ceil_to_nearest(float(best_design_var['bf']), 10)
            self.total_depth = ceil_to_nearest(float(best_design_var['D']), 25)
        else:
            opt_tf_bot = float(best_design_var['tf_bot'])
            self.bottom_flange_thickness = next((t for t in self.bottom_flange_thickness_list if float(t) >= opt_tf_bot), self.bottom_flange_thickness_list[-1])
            
            opt_tf_top = float(best_design_var['tf_top'])
            self.top_flange_thickness = next((t for t in self.top_flange_thickness_list if float(t) >= opt_tf_top), self.top_flange_thickness_list[-1])

            opt_tw = float(best_design_var['tw'])
            self.web_thickness = next((t for t in self.web_thickness_list if float(t) >= opt_tw), self.web_thickness_list[-1])

            self.bottom_flange_width = ceil_to_nearest(float(best_design_var['bf_bot']), 10)
            self.top_flange_width = ceil_to_nearest(float(best_design_var['bf_top']), 10)
            self.total_depth = ceil_to_nearest(float(best_design_var['D']), 25)

        if not is_thick_web:
            opt_t_stiff = float(best_design_var['t_stiff'])
            self.IntStiffThickness = next((t for t in self.int_thickness_list if float(t) >= opt_t_stiff), self.int_thickness_list[-1])
            self.c = ceil_to_nearest(float(best_design_var['c']), 10)
        
        logger.info(f"Optimized values rounded up: D={self.total_depth}, tw={self.web_thickness}, "
                    f"Bf_top={self.top_flange_width}, tf_top={self.top_flange_thickness}, "
                    f"Bf_bot={self.bottom_flange_width}, tf_bot={self.bottom_flange_thickness}, "
                    f"c={self.c}, t_stiff={self.IntStiffThickness}")

        # Run final design check with rounded values and full logging
        self.design_check(design_dictionary, silent=False)
        
        # Clear optimization context
        self._optimization_variable_list = []
        self._optimization_design_dictionary = {}


    def final_format(self, design_dictionary):
        """
        Populates all 'self' attributes needed for the output dock and reports
        after a design (pass or fail) is complete.
        """
        self.result_designation = (f"{int(self.total_depth)} x {int(self.web_thickness)} x "
                                   f"{int(self.bottom_flange_width)} x {int(self.bottom_flange_thickness)} x "
                                   f"{int(self.top_flange_width)} x {int(self.top_flange_thickness)}")
        
        self.result_UR = max(self.moment_ratio or 0, self.shear_ratio or 0, self.deflection_ratio or 0)
        self.section_classification_val = self.section_class
        self.betab = round(self.beta_b_lt, 2) if self.beta_b_lt else 0
        
        if self.total_depth > self.top_flange_thickness + self.bottom_flange_thickness:
             self.effectivearea = Unsymmetrical_I_Section_Properties.calc_area(self, self.total_depth, self.top_flange_width, self.bottom_flange_width, self.web_thickness, self.top_flange_thickness, self.bottom_flange_thickness) / 100
        else:
             self.effectivearea = 0
             
        self.design_moment = round(self.Md / 1e6, 2) if self.Md else 0
        
        if self.support_type == 'Major Laterally Unsupported':
            self.critical_moment = round(self.M_cr / 1e6, 2) if self.M_cr else 0
            self.torsion_cnst = round(self.It / 1e4, 2) if self.It else 0
            self.warping_cnst = round(self.Iw / 1e6, 2) if self.Iw else 0
        
        self.intstiffener_thk = self.IntStiffThickness
        self.longstiffener_thk = self.LongStiffThickness
        self.longstiffener_no = 0
        if self.long_Stiffner == 'Yes and 1 stiffener': self.longstiffener_no = 1
        elif self.long_Stiffner == 'Yes and 2 stiffeners': self.longstiffener_no = 2
        
        self.intstiffener_spacing = self.c
        self.end_panel_stiffener_thickness = self.end_stiffthickness
        
        self.atop, self.abot = self.design_welds_with_strength_web_to_flange(self, self.load.shear_force, self.top_flange_width, self.top_flange_thickness, self.bottom_flange_width, self.bottom_flange_thickness, self.eff_depth, [self.material.fu])
        
        if self.V_d is not None and self.end_stiffthickness > 0:
            self.weld_stiff = self.weld_for_end_stiffener(self, self.end_stiffthickness, self.end_stiffwidth, self.load.shear_force, self.V_d, self.total_depth, self.top_flange_thickness, self.bottom_flange_thickness, self.web_thickness, [self.material.fu])
        else:
            self.weld_stiff = 0 
            
        # self.design_status is set in the calling function (design_check or optimized_method)

    def save_design(self, popup_summary):
        """
        Saves the design report.
        """
        logger.info(" :=========Start Of Design Saving===========")
        # This function would contain the logic to generate and save
        # the design report using the 'self' attributes populated by final_format.
        
    def get_3d_components(self):
        """
        Returns the components for 3D visualization.
        """
        components = []
        # TODO: Implement 3D model generation
        # t3 = ('Model', self.call_3DModel)
        # components.append(t3)
        return components