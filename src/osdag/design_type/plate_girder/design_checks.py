import math
import numpy as np
import logging
from ...Common import *
from ...utils.common import is800_2007
from ...utils.common.common_calculation import round_up
from ...utils.common.Unsymmetrical_Section_Properties import Unsymmetrical_I_Section_Properties
from ...utils.common.component import *

# Set up a logger for this module
logger = logging.getLogger('Osdag.PlateGirderChecks')

class PlateGirderLogic:
    """
    Contains all the core engineering calculations for welded plate girders,
    decoupled from the Osdag UI framework.
    
    This class is intended to be used as a mixin or base class for the
    main PlateGirderWelded class.
    """

    # --- Section & Material Checks ---

    def section_classification(self, design_dictionary):
        self.design_status = False

        flange_class_top = IS800_2007.Table2_i(((self.top_flange_width / 2)), self.top_flange_thickness, self.material.fy, 'Welded')[0]
        flange_class_bottom = IS800_2007.Table2_i(((self.bottom_flange_width / 2)), self.bottom_flange_thickness, self.material.fy, 'Welded')[0]
        web_class = IS800_2007.Table2_iii((self.total_depth - self.top_flange_thickness - self.bottom_flange_thickness), self.web_thickness, self.material.fy)

        if flange_class_bottom == "Slender" or web_class == "Slender" or flange_class_top == 'Slender':
            self.section_class = "Slender"
        else:
            if flange_class_top == KEY_Plastic:
                if web_class == KEY_Plastic:
                    if flange_class_bottom == KEY_Plastic:
                        self.section_class = KEY_Plastic
                    elif flange_class_bottom == KEY_Compact:
                        self.section_class = KEY_Compact
                    else:  # SemiCompact
                        self.section_class = KEY_SemiCompact
                elif web_class == KEY_Compact:
                    if flange_class_bottom in [KEY_Plastic, KEY_Compact]:
                        self.section_class = KEY_Compact
                    else:  # SemiCompact
                        self.section_class = KEY_SemiCompact
                else:  # web SemiCompact
                    self.section_class = KEY_SemiCompact
            elif flange_class_top == KEY_Compact:
                if web_class == KEY_Plastic:
                    if flange_class_bottom in [KEY_Plastic, KEY_Compact]:
                        self.section_class = KEY_Compact
                    else:  # SemiCompact
                        self.section_class = KEY_SemiCompact
                elif web_class == KEY_Compact:
                    if flange_class_bottom in [KEY_Plastic, KEY_Compact]:
                        self.section_class = KEY_Compact
                    else:  # SemiCompact
                        self.section_class = KEY_SemiCompact
                else:  # web SemiCompact
                    self.section_class = KEY_SemiCompact
            else:  # flange_class_top == SemiCompact
                self.section_class = KEY_SemiCompact

        self.Zp_req = self.load.moment * self.gamma_m0 / self.material.fy
        self.effective_length_beam(self, design_dictionary, self.length)

        if self.section_class == 'Slender' and self.web_philosophy == 'Thick Web without ITS':
            return False
        else:
            return True

    def beta_value(self, design_dictionary, section_class):
        self.plast_sec_mod_z = Unsymmetrical_I_Section_Properties.calc_PlasticModulusZ(self, self.total_depth, self.top_flange_width, self.bottom_flange_width,
                                                                                       self.web_thickness, self.top_flange_thickness, self.bottom_flange_thickness, self.epsilon)
        self.elast_sec_mod_z = Unsymmetrical_I_Section_Properties.calc_ElasticModulusZz(self, self.total_depth, self.top_flange_width, self.bottom_flange_width,
                                                                                       self.web_thickness, self.top_flange_thickness, self.bottom_flange_thickness)
        self.Zp_req = self.load.moment * self.gamma_m0 / self.material.fy
        if self.plast_sec_mod_z < self.Zp_req:
            pass # Fails check, but beta value is still calculated
            
        if section_class == KEY_Plastic or section_class == KEY_Compact:
            self.beta_b_lt = 1.0
        else:
            if self.plast_sec_mod_z == 0:
                self.beta_b_lt = 1.0 # Avoid division by zero, set to semi-compact limit
            else:
                self.beta_b_lt = (self.elast_sec_mod_z / self.plast_sec_mod_z)

    def min_web_thickness_thick_web(self, d, tw, eps, stiffener_type, c, silent=False):
        if IS800_2007.cl_8_6_1_1_and_8_6_1_2_web_thickness_check(d, tw, eps, stiffener_type, c):
            return True
        else:
            if not silent:
                logger.error("Web thickness is insufficient.")
            return False

    # --- Moment Capacity Checks ---

    def moment_capacity_laterally_supported(self, V, Zp, Ze, Fy, gamma_m0, D, tw, tf_top, tf_bot, section_class):
        A_vg = (D - tf_top - tf_bot) * tw
        self.V_d = ((A_vg * Fy) / (math.sqrt(3) * gamma_m0))
        if V > 0.6 * self.V_d:  # high shear
            self.Md = self.calc_Mdv(self, V, self.V_d, Zp, Ze, Fy, gamma_m0, D, tw, tf_top, tf_bot)
        else:  # low shear
            self.Md = IS800_2007.cl_8_2_1_2_design_bending_strength(section_class, Zp, Ze, Fy, gamma_m0, self.support_condition)
        
        if self.Md == 0:
            self.moment_ratio = 2.0 # Fail
            return False
            
        self.moment_ratio = self.load.moment / self.Md
        return self.Md >= self.load.moment

    def moment_capacity_laterally_unsupported(self, E, LLT, D,
                                              tf_top, tf_bot, Bf_top, Bf_bot, tw,
                                              LoadingCase, gamma_m0, Fy, shear_force, silent=False):
        if Bf_top == Bf_bot and tf_top == tf_bot:
            yj = 0
        else:
            h = D - (tf_top + tf_bot)
            Ift = (Bf_top * tf_top ** 3) / 12
            Ifc = (Bf_bot * tf_bot ** 3) / 12
            if (Ifc + Ift) == 0:
                beta_f = 0
            else:
                beta_f = Ifc / (Ifc + Ift)
            alpha = 0.8 if beta_f > 0.5 else 1.0
            yj = alpha * (2 * beta_f - 1) * h / 2

        G = 0.769 * 10 ** 5
        Kw = self.get_K_from_warping_restraint(self, self.warping, silent)
        Iy = Unsymmetrical_I_Section_Properties.calc_MomentOfAreaY(self, self.total_depth, self.top_flange_width,
                                                                   self.bottom_flange_width, self.web_thickness,
                                                                   self.top_flange_thickness,
                                                                   self.bottom_flange_thickness)
        self.It = Unsymmetrical_I_Section_Properties.calc_TorsionConstantIt(self, self.total_depth, self.top_flange_width,
                                                                       self.bottom_flange_width, self.web_thickness,
                                                                       self.top_flange_thickness,
                                                                       self.bottom_flange_thickness)
        self.Iw = Unsymmetrical_I_Section_Properties.calc_WarpingConstantIw(self, self.total_depth, self.top_flange_width,
                                                                       self.bottom_flange_width, self.web_thickness,
                                                                       self.top_flange_thickness,
                                                                       self.bottom_flange_thickness)

        # Mcr calc
        yg = D / 2
        yj = self.calc_yj(self, Bf_top, tf_top, Bf_bot, tf_bot, D)
        K_value = 0
        
        if LoadingCase == KEY_DISP_UDL_PIN_PIN_PG:
            K_value = 1.0
            c1, c2, c3 = 1.132, 0.459, 0.525
        elif LoadingCase == KEY_DISP_UDL_FIX_FIX_PG:
            K_value = 0.5
            c1, c2, c3 = 0.712, 0.652, 1.070
        elif LoadingCase == KEY_DISP_PL_PIN_PIN_PG:
            K_value = 1.0
            c1, c2, c3 = 1.365, 0.553, 1.780
        elif LoadingCase == KEY_DISP_PL_FIX_FIX_PG:
            K_value = 0.5
            c1, c2, c3 = 0.938, 0.715, 4.800
        else:
            # Default to UDL PIN_PIN as a fallback
            K_value = 1.0
            c1, c2, c3 = 1.132, 0.459, 0.525
            if not silent:
                logger.warning(f"Invalid Loading Case '{LoadingCase}'. Defaulting to UDL Pin-Pin.")

        if Iy == 0 or E == 0: # Avoid division by zero
             self.M_cr = 0
        elif Bf_top == Bf_bot and tf_top == tf_bot: # Symmetric
            term1 = (math.pi ** 2 * E * Iy) / (LLT ** 2)
            term2 = (self.Iw / Iy)
            term3 = (G * self.It * LLT ** 2) / (math.pi ** 2 * E * Iy)
            self.M_cr = term1 * math.sqrt(term2 + term3)
        else: # Unsymmetric
            term1 = (math.pi ** 2 * E * Iy) / (LLT ** 2)
            bracket_term1 = (K_value / Kw) ** 2 * (self.Iw / Iy)
            bracket_term2 = (G * self.It * LLT ** 2) / (math.pi ** 2 * E * Iy)
            bracket_term3 = (c2 * yg - c3 * yj) ** 2
            bracket = bracket_term1 + bracket_term2 + bracket_term3
            self.M_cr = c1 * term1 * math.sqrt(bracket) - term1 * (c2 * yg - c3 * yj)

        A_vg = (D - tf_top - tf_bot) * tw
        self.V_d = ((A_vg * Fy) / (math.sqrt(3) * gamma_m0))
        
        if self.M_cr <= 0: # Handle non-physical M_cr
             self.Md = 0
        else:
            self.lambda_lt = IS800_2007.cl_8_2_2_1_elastic_buckling_moment(self.beta_b_lt, self.plast_sec_mod_z, self.elast_sec_mod_z, Fy,
                                                                        self.M_cr)

            self.phi_lt = IS800_2007.cl_8_2_2_Unsupported_beam_bending_phi_lt(self.alpha_lt, self.lambda_lt)
            self.X_lt = IS800_2007.cl_8_2_2_Unsupported_beam_bending_stress_reduction_factor(self.phi_lt, self.lambda_lt)
            self.fbd_lt = IS800_2007.cl_8_2_2_Unsupported_beam_bending_compressive_stress(self.X_lt, Fy, self.gamma_m0)
            self.Md = IS800_2007.cl_8_2_2_Unsupported_beam_bending_strength(self.plast_sec_mod_z, self.elast_sec_mod_z, self.fbd_lt,
                                                                            self.section_class)
        
        if shear_force > 0.6 * self.V_d:  # high shear
            self.Md = self.calc_Mdv_lat_unsupported(self, self.load.shear_force, self.V_d, self.plast_sec_mod_z,
                                                    self.elast_sec_mod_z, self.material.fy, self.gamma_m0,
                                                    self.total_depth, self.web_thickness, self.top_flange_thickness,
                                                    self.bottom_flange_thickness, self.Md)
        # Low shear check is implicit, self.Md is already set
        
        if self.Md == 0:
            self.moment_ratio = 2.0 # Fail
            return False
            
        self.moment_ratio = self.load.moment / self.Md
        return self.Md >= self.load.moment

    def calc_Mdv(self, V, Vd, Zp, Ze, Fy, gamma_m0, D, tw, tf_top, tf_bot):
        beta = (2 * V / Vd - 1) ** 2
        d = D - (tf_top + tf_bot)
        Aw = d * tw
        Zfd = Zp - (Aw * D / 4)
        Mfd = Zfd * Fy / gamma_m0
        Md = Zp * Fy / gamma_m0
        Mdv = Md - beta * (Md - Mfd)
        Mdv_limit = (1.2 * Ze * Fy) / gamma_m0
        return round(min(Mdv, Mdv_limit), 2)

    def calc_Mdv_lat_unsupported(self, V, Vd, Zp, Ze, Fy, gamma_m0, D, tw, tf_top, tf_bot, Md_low_shear):
        beta = (2 * V / Vd - 1) ** 2
        d = D - (tf_top + tf_bot)
        Aw = d * tw
        Zfd = Zp - (Aw * D / 4)
        Mfd = Zfd * Fy / gamma_m0
        Mdv = Md_low_shear - beta * (Md_low_shear - Mfd) # Md_low_shear is the M_d calculated using 8.2.2
        Mdv_limit = (1.2 * Ze * Fy) / gamma_m0
        return round(min(Mdv, Mdv_limit), 2)

    # --- Shear Capacity Checks ---

    def shear_capacity_laterally_supported_thick_web(self, Fy, gamma_m0, D, tw, tf_top, tf_bot):
        A_vg = (D - tf_top - tf_bot) * tw
        if A_vg <= 0:
            self.V_d = 0
            self.shear_ratio = 2.0 # Fail
            return False
            
        self.V_d = ((A_vg * Fy) / (math.sqrt(3) * gamma_m0))
        self.shear_ratio = self.load.shear_force / self.V_d
        return self.V_d >= self.load.shear_force

    def shear_buckling_check_simple_postcritical(self, eff_depth, D, tf_top, tf_bot, tw, V, c=0):
        if eff_depth <= 0: # Invalid geometry
            self.V_cr = 0
            self.shear_ratio = 2.0
            return False
            
        A_vg = eff_depth * tw
        if c == 0 or (c / eff_depth) >= 1:
            K_v = 5.35 + 4 / ((c / eff_depth) ** 2) if c != 0 else 5.35
        else:
            K_v = 4 + 5.35 / (c / eff_depth) ** 2
        
        E = self.material.modulus_of_elasticity
        mu = 0.3
        tau_crc = IS800_2007.cl_8_4_2_2_tau_crc_Simple_postcritical(K_v, E, mu, eff_depth, self.web_thickness)
        lambda_w = IS800_2007.cl_8_4_2_2_lambda_w_Simple_postcritical(self.material.fy, tau_crc)
        tau_b = IS800_2007.cl_8_4_2_2_tau_b_Simple_postcritical(lambda_w, self.material.fy)
        self.V_cr = IS800_2007.cl_8_4_2_2_Vcr_Simple_postcritical(tau_b, A_vg)
        
        if self.V_cr > V:
            self.shear_ratio = max(self.load.shear_force / self.V_cr, self.shear_ratio)
            return True
        else:
            self.shear_ratio = max(self.load.shear_force / self.V_cr, self.shear_ratio)
            return False # Fails simple check, will require end stiffener check

    def shear_buckling_check_tension_field(self, eff_depth, D, tf_top, tf_bot, tw, c=0):
        if eff_depth <= 0: # Invalid geometry
            self.V_tf = 0
            self.shear_ratio = 2.0
            return False
            
        A_vg = (D - tf_top - tf_bot) * tw
        if c == 0 or (c / eff_depth) >= 1:
            K_v = 5.35 + 4 / ((c / eff_depth) ** 2) if c != 0 else 5.35
        else:
            K_v = 4 + 5.35 / (c / eff_depth) ** 2
        
        E = self.material.modulus_of_elasticity
        mu = 0.3
        tau_crc = IS800_2007.cl_8_4_2_2_tau_crc_Simple_postcritical(K_v, E, mu, eff_depth, self.web_thickness)
        lambda_w = IS800_2007.cl_8_4_2_2_lambda_w_Simple_postcritical(self.material.fy, tau_crc)
        tau_b = IS800_2007.cl_8_4_2_2_tau_b_Simple_postcritical(lambda_w, self.material.fy)
        self.V_cr = IS800_2007.cl_8_4_2_2_Vcr_Simple_postcritical(tau_b, A_vg)
        
        Nf = self.load.moment / (eff_depth + (tf_top + tf_bot) / 2)
        _phi, _M_fr_t, _M_fr_b, _s_t, _s_b, _w_tf, _sai, _fv, self.V_tf = IS800_2007.cl_8_4_2_2_TensionField_unequal_Isection(
            c, eff_depth, self.web_thickness, self.material.fy, self.top_flange_width,
            self.top_flange_thickness, self.bottom_flange_width, self.bottom_flange_thickness,
            Nf, self.gamma_m0, A_vg, tau_b
        )
        
        self.shear_ratio = max(self.load.shear_force / self.V_tf, self.shear_ratio)
        return self.V_tf >= self.load.shear_force

    # --- Web Buckling / Crippling Checks ---

    def web_buckling_laterally_supported_thick_web(self, Fy, gamma_m0, D, tw, tf_top, tf_bot, E, b1):
        self.eff_depth = D - (tf_bot + tf_top)
        if self.eff_depth <= 0: # Invalid geometry
             self.web_buckling_ratio = 2.0
             self.shear_ratio = max(self.web_buckling_ratio, self.shear_ratio)
             return False
             
        n1 = self.eff_depth / 2 # Dispersion
        Ac = (b1 + n1) * tw
        slenderness_input = 2.5 * self.eff_depth / tw
        self.fcd = IS800_2007.cl_7_1_2_1_design_compressisive_stress_plategirder(Fy, gamma_m0, slenderness_input, E)
        Critical_buckling_load = round(Ac * self.fcd, 2)
        
        if Critical_buckling_load == 0:
            self.web_buckling_ratio = 2.0 # Fail
        else:
            self.web_buckling_ratio = self.load.shear_force / Critical_buckling_load
            
        self.shear_ratio = max(self.web_buckling_ratio, self.shear_ratio)
        return Critical_buckling_load >= self.load.shear_force

    def web_crippling_laterally_supported_thick_web(self, Fy, gamma_m0, tw, tf_top, b1, silent=False):
        try:
            web_height = self.total_depth - self.top_flange_thickness - self.bottom_flange_thickness
            if web_height <= 0: # Invalid geometry
                return False
            return self.check_web_crippling(b1, tw, Fy, web_height, silent)
        except Exception as e:
            if not silent:
                logger.error(f"Error in web crippling check: {str(e)}")
            return False

    def check_web_crippling(self, N, tw, fy, d, silent=False):
        """
        Check web crippling as per IS 800:2007 Section 8.7.2.
        N = b1 (bearing length)
        d = clear depth of web
        """
        try:
            if any(val <= 0 for val in [N, tw, fy, d]):
                if not silent:
                    # Use the class-level logger if available, else print
                    log_msg = "Invalid input parameters (<= 0) for web crippling check"
                    try:
                        logger.warning(log_msg)
                    except NameError:
                        print(f"WARNING: {log_msg}")
                return False

            k1 = 3.25  # For end reactions (Clause 8.7.2.1)
            k2 = 0.15  # For end reactions (Clause 8.7.2.1)
            E = self.material.modulus_of_elasticity

            P_w = (k1 * k2 * tw * tw * math.sqrt(fy * E)) * (1 + (N / d))
            P_w = P_w / self.gamma_m0

            if d / tw > 200:
                if not silent:
                    log_msg = "Web slenderness ratio (d/tw) exceeds 200. Additional stiffening may be required."
                    try:
                        logger.warning(log_msg)
                    except NameError:
                        print(f"WARNING: {log_msg}")

            if P_w >= self.load.shear_force:
                return True
            else:
                if not silent:
                    log_msg = f"Web crippling resistance ({P_w:.2f} N) is less than factored load ({self.load.shear_force:.2f} N)"
                    try:
                        logger.warning(log_msg)
                    except NameError:
                        print(f"WARNING: {log_msg}")
                return False

        except Exception as e:
            if not silent:
                log_msg = f"Error in web crippling calculation: {str(e)}"
                try:
                    logger.error(log_msg)
                except NameError:
                    print(f"ERROR: {log_msg}")
            return False

    # --- Stiffener Design ---

    def end_panel_stiffener_calc(self,
                                 Bf_top, Bf_bot, tw, fy, gamma_m0, d,
                                 tf_top, total_depth, effective_length, tf_bot, E, eps, c, silent=False
                                 ):
        if d <= 0: # Invalid geometry
            self.endshear_ratio = 2.0
            self.shear_ratio = max(self.endshear_ratio, self.shear_ratio)
            return False
            
        A_vg = d * tw
        if c is None or c == 0 or c == 'NA':
            c = d  # Default to web depth as recommended

        c = float(c) # Ensure c is float for calculations
        
        if (c / d) < 1:
            K_v = 4 + 5.35 / (c / d) ** 2
        else:
            K_v = 5.35 + 4 / (c / d) ** 2
        
        mu = 0.3
        tau_crc = IS800_2007.cl_8_4_2_2_tau_crc_Simple_postcritical(K_v, E, mu, d, tw)
        lambda_w = IS800_2007.cl_8_4_2_2_lambda_w_Simple_postcritical(fy, tau_crc)
        tau_b = IS800_2007.cl_8_4_2_2_tau_b_Simple_postcritical(lambda_w, fy)
        self.V_cr = IS800_2007.cl_8_4_2_2_Vcr_Simple_postcritical(tau_b, A_vg)
        
        Nf = self.load.moment / d # Approx
        
        # This part of the stiffener check (tension-field-related forces)
        # seems misplaced for a general end-panel check per 8.7.4
        # Re-focusing on 8.7.4: Check buckling (8.7.4.1) and bearing (8.7.4.2)
        
        thickness_list = ['8', '10', '12', '14', '16', '18', '20', '22', '25', '28', '32', '36', '40', '45', '50', '56', '63', '75', '80', '90', '100', '110', '120']
        
        if not self.int_thickness_list:
            if not silent:
                logger.error("Intermediate stiffener thickness list is empty.")
            return False

        found_stiffener = False
        for stiff_thick_str in thickness_list:
            self.end_stiffthickness = float(stiff_thick_str)
            
            # Use the instance's end_stiffwidth
            current_stiff_width = self.end_stiffwidth
            if current_stiff_width <= 0: # Ensure valid width
                current_stiff_width = (min(Bf_top, Bf_bot) - tw) / 2
                if current_stiff_width <= 0: current_stiff_width = 50 # fallback
            
            max_outstand = 14 * self.end_stiffthickness * eps # 8.7.1.2
            if current_stiff_width > max_outstand:
                current_stiff_width = max_outstand
            
            min_thickness_check = current_stiff_width / 16 # 8.7.1.3
            if self.end_stiffthickness < min_thickness_check:
                continue
                
            # 8.7.4.2: Effective section for bearing
            web_contrib_length_bearing = min(20 * tw, d / 2)
            A_q_bearing = (2 * current_stiff_width * self.end_stiffthickness) + (web_contrib_length_bearing * tw)
            
            # 8.7.4.1: Effective section for buckling
            web_contrib_length_buckling = min(20 * tw, c/2) # Use 'c' for buckling
            A_q_buckling = (2 * current_stiff_width * self.end_stiffthickness) + (web_contrib_length_buckling * tw)
            
            # I_x for buckling (8.7.4.1)
            # I of stiffener plates about web centerline
            I_x_stiff = 2 * ( (self.end_stiffthickness * (current_stiff_width**3)) / 12 + (current_stiff_width * self.end_stiffthickness) * (tw/2 + current_stiff_width/2)**2 )
            I_x_web = (web_contrib_length_buckling * (tw**3)) / 12
            I_x = I_x_stiff + I_x_web
            
            if A_q_buckling == 0: continue # Avoid division by zero
            r_x = math.sqrt(I_x / A_q_buckling)
            Le = 0.7 * d # 8.7.4.1
            KL_r = Le / r_x if r_x > 0 else float('inf')
            
            fcd = IS800_2007.cl_7_1_2_1_design_compressisive_stress_plategirder(fy, gamma_m0, KL_r, E)
            Pd = A_q_buckling * fcd # Buckling Resistance
            self.Critical_buckling_resistance = Pd

            # Bearing Capacity (8.7.4.2)
            P_bf = A_q_bearing * fy / gamma_m0
            
            # The force on the stiffener is the shear force V
            force_on_stiffener = self.load.shear_force
            
            if Pd == 0 or P_bf == 0:
                self.endshear_ratio = 2.0 # Fail
                continue

            # Check 1: Buckling Resistance (8.7.4.1)
            buckling_ratio = force_on_stiffener / Pd
            # Check 2: Bearing Capacity (8.7.4.2)
            bearing_ratio = force_on_stiffener / P_bf
            
            self.endshear_ratio = max(buckling_ratio, bearing_ratio)
            
            if self.endshear_ratio <= 1.0:
                found_stiffener = True
                break
            else:
                continue
                
        self.shear_ratio = max(self.endshear_ratio, self.shear_ratio)
        
        if not found_stiffener:
            self.end_stiffthickness = 0 # Reset if check fails
            return False
        
        return True

    def shear_buckling_check_intermediate_stiffener(
        self, d, tw, c, e, IntStiffThickness, IntStiffenerWidth, V_ed, gamma_m0, fy, E
    ):
        if d <= 0 or c <= 0: return False # Invalid geometry
        
        A_vg = d * tw
        c = float(c) # Ensure c is float
        
        if (c / d) < 1:
            K_v = 4 + 5.35 / (c / d) ** 2
        else:
            K_v = 5.35 + 4 / (c / d) ** 2
        
        mu = 0.3
        tau_crc = IS800_2007.cl_8_4_2_2_tau_crc_Simple_postcritical(K_v, E, mu, d, tw)
        lambda_w = IS800_2007.cl_8_4_2_2_lambda_w_Simple_postcritical(fy, tau_crc)
        tau_b = IS800_2007.cl_8_4_2_2_tau_b_Simple_postcritical(lambda_w, fy)
        self.V_cr = IS800_2007.cl_8_4_2_2_Vcr_Simple_postcritical(tau_b, A_vg)
        
        # Check stiffener rigidity (8.7.2.2)
        cd_ratio = c / d
        I_min_global = 0.75 * d * tw**3 if cd_ratio >= math.sqrt(2) else (1.5 * d**3 * tw**3) / (c**2)

        max_outstand = 14 * IntStiffThickness * e # 8.7.1.2
        if IntStiffenerWidth > max_outstand:
            IntStiffenerWidth = max_outstand

        # I of stiffener plates about web centerline
        I_s_plates = 2 * ( (IntStiffThickness * (IntStiffenerWidth**3)) / 12 + (IntStiffenerWidth * IntStiffThickness) * (tw/2 + IntStiffenerWidth/2)**2 )
        I_s = I_s_plates
        
        if I_s < I_min_global:
             return False # Failed global inertia check

        # Check stiffener buckling (8.7.2.4)
        F_q = (V_ed - self.V_cr) # This is the force to be resisted
        if F_q <= 0: 
             self.Critical_buckling_resistance = float('inf')
             return True # Passes axial check, no force

        A_s = 2 * IntStiffenerWidth * IntStiffThickness
        web_contrib_length_buckling = min(20 * tw, c/2)
        A_x = A_s + (web_contrib_length_buckling * tw) # Effective area for buckling
        
        I_x_web = (web_contrib_length_buckling * (tw**3)) / 12
        I_x = I_s + I_x_web # Total I about web centerline
        
        if A_x == 0: return False # Avoid division by zero
        
        r_x = math.sqrt(I_x / A_x)
        Le = d # 8.7.2.4
        slenderness_input = Le / r_x if r_x > 0 else float('inf')
        
        fcd = IS800_2007.cl_7_1_2_1_design_compressisive_stress_plategirder(
            fy, gamma_m0, slenderness_input, E
        )
        Pd = round(A_x * fcd, 2)
        
        self.shear_ratio = max(F_q / Pd if Pd > 0 else 2.0, self.shear_ratio)
        self.Critical_buckling_resistance = Pd

        return F_q <= Pd


    def tension_field_end_stiffener(self, *args, **kwargs):
        # The logic for an end stiffener is the same regardless of
        # simple post-critical or tension field method.
        # IS 8.7.4 applies.
        return self.end_panel_stiffener_calc(*args, **kwargs)


    def tension_field_intermediate_stiffener(self, *args, **kwargs):
        # IS 8.7.2.4 applies for intermediate stiffeners in tension-field design
        # This logic is captured in shear_buckling_check_intermediate_stiffener
        return self.shear_buckling_check_intermediate_stiffener(*args, **kwargs)

    def design_longitudinal_stiffeners(self, d, tw, c, eps_w, second_stiffener=False, silent=False):
        if d <= 0 or c <= 0 or tw <= 0: return False # Invalid geometry
        
        c = float(c)
        tw = float(tw)
        d_na = Unsymmetrical_I_Section_Properties.calc_centroid(self, self.total_depth, self.top_flange_width, self.bottom_flange_width, self.web_thickness, self.top_flange_thickness, self.bottom_flange_thickness)

        self.x1 = int(round(d_na * 0.2, 0)) # 0.2 * d_c (dist from comp flange)
        self.x2 = int(round(d_na, 0)) # at NA

        # 8.7.2.5 (b) - Rigidity of stiffener at 0.2d_c
        I1_min = 4.0 * d * tw ** 3 
        
        # 8.7.3.2 - Rigidity of stiffener at NA
        I2_min = d * tw ** 3 # Note: IS 800:2007 8.7.3.2
        
        # Provided I_s (Moment of inertia of the stiffener about the web)
        width_ls = self.eff_width_longitudnal
        thick_ls = self.LongStiffThickness
        if width_ls <= 0 or thick_ls <= 0: return False # No stiffener
        
        A_ls = width_ls * thick_ls
        # I of stiffener about its own centroid + A*d^2 (Parallel Axis Theorem)
        # Assuming stiffener is centered on web (e.g. two plates)
        # This check is complex and needs clear definition of stiffener geometry.
        # Using a simplification: I of a plate about its edge.
        Is_provided_1 = (thick_ls * (width_ls**3)) / 3 # Simplified
        Is_provided_2 = (thick_ls * (width_ls**3)) / 3 # Simplified

        if second_stiffener is False:
            return Is_provided_1 >= I1_min
        else:
            return Is_provided_1 >= I1_min and Is_provided_2 >= I2_min

    # --- Serviceability (Deflection) ---

    def deflection_from_moment_kNm_mm(self, M_kNm, L, E, I, case, silent=False):
        M = M_kNm
        if E == 0 or I == 0:
            return float('inf') # Fail
        
        pref = M * L ** 2 / (E * I * 1.5)
        
        if case == KEY_DISP_UDL_PIN_PIN_PG:
            return (5 / 48) * pref
        elif case == KEY_DISP_UDL_FIX_FIX_PG:
            return (1 / 32) * pref
        elif case == KEY_DISP_PL_PIN_PIN_PG:
            return (1 / 12) * pref
        elif case == KEY_DISP_PL_FIX_FIX_PG:
            return (1 / 24) * pref
        else:
            if not silent:
                logger.warning(f"Unknown deflection case '{case}'. Defaulting to 'simple_udl'.")
            return (5 / 48) * pref

    def evaluate_deflection_kNm_mm(self, M_kNm, L, E, case, criteria, silent=False):
        I = Unsymmetrical_I_Section_Properties.calc_MomentOfAreaZ(self, self.total_depth, self.top_flange_width, self.bottom_flange_width, self.web_thickness, self.top_flange_thickness, self.bottom_flange_thickness)
        
        if I == 0: # Avoid division by zero
             self.deflection_ratio = 2.0 # Fail
             return False
             
        delta = self.deflection_from_moment_kNm_mm(self, M_kNm, L, E, I, case, silent)

        if criteria == 'NA':
            if not silent:
                logger.info("Deflection criteria is 'NA'. Skipping deflection check.")
            self.deflection_ratio = 0
            return True # Pass if no criteria

        try:
            n = float(criteria)
            if n <= 0:
                 if not silent:
                     logger.error(f"Invalid deflection criteria 'L/{n}'. Skipping check.")
                 self.deflection_ratio = 0
                 return True
        except ValueError:
            if not silent:
                logger.error(f"Could not parse deflection criteria 'L/{criteria}'. Skipping check.")
            self.deflection_ratio = 0
            return True

        allowable = L / n
        if allowable == 0:
            self.deflection_ratio = 2.0 # Fail
            return False
            
        ok = (delta <= allowable)
        self.deflection_ratio = delta / allowable
        return ok

    # --- Weld Design ---

    def shear_stress_unsym_I(self, V_ed, b_ft, t_ft, b_fb, t_fb, t_w, h_w):
        if h_w <= 0: return {'q_top_kN_per_mm': 0, 'q_bot_kN_per_mm': 0}
        
        A_t = b_ft * t_ft
        A_b = b_fb * t_fb
        A_w = t_w * h_w
        A = A_t + A_b + A_w
        if A == 0: return {'q_top_kN_per_mm': 0, 'q_bot_kN_per_mm': 0}

        y_b = t_fb / 2
        y_w = t_fb + h_w / 2
        y_t = t_fb + h_w + t_ft / 2

        y_na = (A_b * y_b + A_w * y_w + A_t * y_t) / A

        I_b = b_fb * t_fb ** 3 / 12 + A_b * (y_b - y_na) ** 2
        I_w = t_w * h_w ** 3 / 12 + A_w * (y_w - y_na) ** 2
        I_t = b_ft * t_ft ** 3 / 12 + A_t * (y_t - y_na) ** 2
        I_z = I_b + I_w + I_t
        if I_z == 0: return {'q_top_kN_per_mm': 0, 'q_bot_kN_per_mm': 0}

        Q_bot = A_b * abs(y_na - y_b)
        Q_top = A_t * abs(y_t - y_na)

        q_bot = V_ed * Q_bot / I_z
        q_top = V_ed * Q_top / I_z

        return {
            'q_top_kN_per_mm': q_top,
            'q_bot_kN_per_mm': q_bot,
        }

    def weld_leg_from_q_with_cl10(self, q_kN_per_mm, ultimate_stresses, fabrication='shop'):
        f_wd = IS800_2007.cl_10_5_7_1_1_fillet_weld_design_stress(ultimate_stresses)
        if f_wd == 0: return 0
        q_N_per_mm = q_kN_per_mm
        t_throat = q_N_per_mm / f_wd
        return t_throat * math.sqrt(2)

    def design_welds_with_strength_web_to_flange(self, V_ed, b_ft, t_ft, b_fb, t_fb, t_w, h_w, ultimate_stresses):
        sf = self.shear_stress_unsym_I(self, V_ed, b_ft, t_ft, b_fb, t_fb, t_w, h_w)
        min_weld_legtop = IS800_2007.cl_10_5_2_3_min_weld_size(t_ft, t_w)
        min_weld_legbot = IS800_2007.cl_10_5_2_3_min_weld_size(t_fb, t_w)
        max_weld_legtop = IS800_2007.cl_10_5_3_1_max_weld_throat_thickness(t_ft, t_w)
        max_weld_legbot = IS800_2007.cl_10_5_3_1_max_weld_throat_thickness(t_fb, t_w)

        a_top_calc = self.weld_leg_from_q_with_cl10(self, sf['q_top_kN_per_mm'], ultimate_stresses)
        a_bot_calc = self.weld_leg_from_q_with_cl10(self, sf['q_bot_kN_per_mm'], ultimate_stresses)

        a_top = round_up(max(a_top_calc, min_weld_legtop), 1)
        a_bot = round_up(max(a_bot_calc, min_weld_legbot), 1)

        # Re-check against max
        a_top = min(a_top, max_weld_legtop)
        a_bot = min(a_bot, max_weld_legbot)

        return a_top, a_bot

    def weld_for_end_stiffener(self, t_st, b_st, V_ed, V_unstf, D, t_ft, t_fb, tw, ultimate_stresses):
        L_weld = D - t_ft - t_fb
        if L_weld <= 0: return 0 # Invalid geometry
        if b_st <= 0: return 0 # Avoid division by zero
        if V_unstf is None: V_unstf = 0 # Ensure V_unstf is a number

        q1 = tw ** 2 / (5 * b_st)
        delta_V = max(V_ed - V_unstf, 0)
        q2 = delta_V / L_weld
        q_tot = q1 + q2
        q_each = q_tot / 2

        min_weld_leg = IS800_2007.cl_10_5_2_3_min_weld_size(t_st, tw)
        max_weld_leg = IS800_2007.cl_10_5_3_1_max_weld_throat_thickness(t_st, tw)

        weld_stiff_calc = self.weld_leg_from_q_with_cl10(self, q_each, ultimate_stresses)
        weld_stiff = round_up(max(weld_stiff_calc, min_weld_leg), 1)
        weld_stiff = min(weld_stiff, max_weld_leg)
        
        return weld_stiff

    # --- Optimization Helper Functions ---

    def generate_first_particle(self, L, M, fy, is_thick_web, is_symmetric, k=67):
        D_empirical = L / 25
        d_opt = ((M * k) / fy) ** (1/3)
        D_final = max(D_empirical, d_opt, 200) # Ensure min depth

        bf_top = bf_bot = bf = 0.3 * D_final
        e = math.sqrt(250 / fy)
        tf_top = max(bf_top / 24, bf_top / (8.4 * e), 6) # Ensure min thick
        tf_bot = max(bf_bot / 24, bf_bot / (8.4 * e), 6) # Ensure min thick
        tf = max(bf / 24, bf_bot / (8.4 * e), 6) # Ensure min thick

        d = D_final - 2 * tf # Approx
        if d <= 0: d = D_final * 0.8 # Fallback
        
        if is_thick_web:
            tw = max(d / 200, d / (84 * e), 8)
        else:
            tw = max(d / 200, d / (105 * e), 8)

        c = d # Initial guess for spacing
        t_stiff = 6
        
        varlst = []
        if is_symmetric:
            if is_thick_web:
                varlst += [tf, tw, bf, D_final]
            else:
                varlst += [tf, tw, bf, D_final, c, t_stiff]
        else:
            if is_thick_web:
                varlst += [tf_top, tf_bot, tw, bf_top, bf_bot, D_final]
            else:
                varlst += [tf_top, tf_bot, tw, bf_top, bf_bot, D_final, c, t_stiff]
        return varlst

    def build_variable_structure(self, is_thick_web=True, is_symmetric=True):
        variables = []
        if is_symmetric:
            variables += ['tf', 'tw', 'bf', 'D']
        else:
            variables += ['tf_top', 'tf_bot', 'tw', 'bf_top', 'bf_bot', 'D']

        if not is_thick_web:
            variables += ['c', 't_stiff']
        return variables

    def get_bounds(self, variable_list):
        bounds_map = {
            'tf': (6, 100), 'tf_top': (6, 100), 'tf_bot': (6, 100),
            'tw': (6, 40),
            'bf': (100, 1000), 'bf_top': (100, 1000), 'bf_bot': (100, 1000),
            'D': (200, 3000), # Increased upper bound for D
            'c': (200, 3000), # c > 0.2*d
            't_stiff': (6, 40)
        }
        lower = [bounds_map[v][0] for v in variable_list]
        upper = [bounds_map[v][1] for v in variable_list]
        return (np.array(lower), np.array(upper))

    def assign_particle_to_section(self, particle, variable_list):
        """Assigns particle values to a temporary Section object AND self."""
        from .pso_solver import Section # Local import
        sec = Section()
        
        for name, value in zip(variable_list, particle):
            setattr(sec, name, value)

        if 'tf' in variable_list:
            sec.tf_top = sec.tf_bot = sec.tf
            sec.bf_top = sec.bf_bot = sec.bf
        
        # Sync self attributes for calculation methods
        self.top_flange_thickness = sec.tf_top
        self.bottom_flange_thickness = sec.tf_bot
        self.web_thickness = sec.tw
        self.top_flange_width = sec.bf_top
        self.bottom_flange_width = sec.bf_bot
        self.total_depth = sec.D
        
        if self.total_depth > self.top_flange_thickness + self.bottom_flange_thickness:
            self.eff_depth = sec.D - sec.tf_top - sec.tf_bot
        else:
            self.eff_depth = 0 # Invalid geometry
        
        if self.web_thickness > 0:
            self.IntStiffnerwidth = min(self.top_flange_width, self.bottom_flange_width) - self.web_thickness / 2 - 10
            self.end_stiffwidth = self.IntStiffnerwidth
        else:
            self.IntStiffnerwidth = 0
            self.end_stiffwidth = 0
        
        if 'c' in variable_list:
            self.c = sec.c
            self.IntStiffThickness = sec.t_stiff
        else:
            self.c = 'NA' # Important for thick web
            self.IntStiffThickness = 0
        
        return sec

    # --- Other Helpers ---
    
    def effective_length_beam(self, design_dictionary, length):
        if design_dictionary[KEY_LENGTH_OVERWRITE] == 'NA':
            self.effective_length = IS800_2007.cl_8_3_1_EffLen_Simply_Supported(
                Torsional=self.torsional_res, Warping=self.warping,
                length=length, depth=(self.total_depth / 1000), load=self.loading_condition
            )
        else:
            try:
                factor = float(design_dictionary[KEY_LENGTH_OVERWRITE])
                if factor <= 0:
                    raise ValueError("Factor must be positive")
                self.effective_length = length * factor
            except (ValueError, TypeError):
                logger.warning("Invalid Effective Length Parameter. Defaulting to 1.0")
                design_dictionary[KEY_LENGTH_OVERWRITE] = '1.0'
                self.effective_length = length

    def get_K_from_warping_restraint(self, warping_condition, silent=False):
        """
        Return effective length factor K based on exact warping restraint description (IS 800:2007, Clause E.1).
        """
        if warping_condition == "Both flanges fully restrained":
            return 0.5
        elif warping_condition == "Compression flange fully restrained":
            return 0.7
        elif warping_condition == "Compression flange partially restrained":
            return 0.85
        elif warping_condition == "Warping not restrained in both flanges":
            return 1.0
        else:
            if not silent:
                logger.warning(f"Invalid warping restraint '{warping_condition}'. Defaulting to 1.0")
            return 1.0

    def calc_yj(self, Bf_top, tf_top, Bf_bot, tf_bot, D):
        if Bf_top == Bf_bot and tf_top == tf_bot:
            return 0  # symmetric section
        h = D - (tf_top + tf_bot)
        if h <= 0: return 0 # Invalid geometry
        
        Ift = (Bf_top * tf_top**3) / 12
        Ifc = (Bf_bot * tf_bot**3) / 12
        if (Ifc + Ift) == 0: return 0 # Avoid division by zero
        
        beta_f = Ifc / (Ifc + Ift)
        alpha = 0.8 if beta_f > 0.5 else 1.0
        yj = alpha * (2 * beta_f - 1) * h / 2
        return yj