%% InDEL analysis(v.0.1)
% By V. Pinheiro

% It starts from the output of the selection routine, and requires an
% Enrichment matrix as an input.


%% Simplified 3mers

TwomerPlus = {'AA_' ;	'AC_' ;	'AD_' ;	'AE_' ;	'AF_' ;	'AG_' ;	'AH_' ;	'AI_' ;	'AK_' ;	'AL_' ;	'AM_' ;	'AN_' ;	'AP_' ;	'AQ_' ;	'AR_' ;	'AS_' ;	'AT_' ;	'AV_' ;	'AW_' ;	'AY_' ;
'CA_' ;	'CC_' ;	'CD_' ;	'CE_' ;	'CF_' ;	'CG_' ;	'CH_' ;	'CI_' ;	'CK_' ;	'CL_' ;	'CM_' ;	'CN_' ;	'CP_' ;	'CQ_' ;	'CR_' ;	'CS_' ;	'CT_' ;	'CV_' ;	'CW_' ;	'CY_' ;
'DA_' ;	'DC_' ;	'DD_' ;	'DE_' ;	'DF_' ;	'DG_' ;	'DH_' ;	'DI_' ;	'DK_' ;	'DL_' ;	'DM_' ;	'DN_' ;	'DP_' ;	'DQ_' ;	'DR_' ;	'DS_' ;	'DT_' ;	'DV_' ;	'DW_' ;	'DY_' ;
'EA_' ;	'EC_' ;	'ED_' ;	'EE_' ;	'EF_' ;	'EG_' ;	'EH_' ;	'EI_' ;	'EK_' ;	'EL_' ;	'EM_' ;	'EN_' ;	'EP_' ;	'EQ_' ;	'ER_' ;	'ES_' ;	'ET_' ;	'EV_' ;	'EW_' ;	'EY_' ;
'FA_' ;	'FC_' ;	'FD_' ;	'FE_' ;	'FF_' ;	'FG_' ;	'FH_' ;	'FI_' ;	'FK_' ;	'FL_' ;	'FM_' ;	'FN_' ;	'FP_' ;	'FQ_' ;	'FR_' ;	'FS_' ;	'FT_' ;	'FV_' ;	'FW_' ;	'FY_' ;
'GA_' ;	'GC_' ;	'GD_' ;	'GE_' ;	'GF_' ;	'GG_' ;	'GH_' ;	'GI_' ;	'GK_' ;	'GL_' ;	'GM_' ;	'GN_' ;	'GP_' ;	'GQ_' ;	'GR_' ;	'GS_' ;	'GT_' ;	'GV_' ;	'GW_' ;	'GY_' ;
'HA_' ;	'HC_' ;	'HD_' ;	'HE_' ;	'HF_' ;	'HG_' ;	'HH_' ;	'HI_' ;	'HK_' ;	'HL_' ;	'HM_' ;	'HN_' ;	'HP_' ;	'HQ_' ;	'HR_' ;	'HS_' ;	'HT_' ;	'HV_' ;	'HW_' ;	'HY_' ;
'IA_' ;	'IC_' ;	'ID_' ;	'IE_' ;	'IF_' ;	'IG_' ;	'IH_' ;	'II_' ;	'IK_' ;	'IL_' ;	'IM_' ;	'IN_' ;	'IP_' ;	'IQ_' ;	'IR_' ;	'IS_' ;	'IT_' ;	'IV_' ;	'IW_' ;	'IY_' ;
'KA_' ;	'KC_' ;	'KD_' ;	'KE_' ;	'KF_' ;	'KG_' ;	'KH_' ;	'KI_' ;	'KK_' ;	'KL_' ;	'KM_' ;	'KN_' ;	'KP_' ;	'KQ_' ;	'KR_' ;	'KS_' ;	'KT_' ;	'KV_' ;	'KW_' ;	'KY_' ;
'LA_' ;	'LC_' ;	'LD_' ;	'LE_' ;	'LF_' ;	'LG_' ;	'LH_' ;	'LI_' ;	'LK_' ;	'LL_' ;	'LM_' ;	'LN_' ;	'LP_' ;	'LQ_' ;	'LR_' ;	'LS_' ;	'LT_' ;	'LV_' ;	'LW_' ;	'LY_' ;
'MA_' ;	'MC_' ;	'MD_' ;	'ME_' ;	'MF_' ;	'MG_' ;	'MH_' ;	'MI_' ;	'MK_' ;	'ML_' ;	'MM_' ;	'MN_' ;	'MP_' ;	'MQ_' ;	'MR_' ;	'MS_' ;	'MT_' ;	'MV_' ;	'MW_' ;	'MY_' ;
'NA_' ;	'NC_' ;	'ND_' ;	'NE_' ;	'NF_' ;	'NG_' ;	'NH_' ;	'NI_' ;	'NK_' ;	'NL_' ;	'NM_' ;	'NN_' ;	'NP_' ;	'NQ_' ;	'NR_' ;	'NS_' ;	'NT_' ;	'NV_' ;	'NW_' ;	'NY_' ;
'PA_' ;	'PC_' ;	'PD_' ;	'PE_' ;	'PF_' ;	'PG_' ;	'PH_' ;	'PI_' ;	'PK_' ;	'PL_' ;	'PM_' ;	'PN_' ;	'PP_' ;	'PQ_' ;	'PR_' ;	'PS_' ;	'PT_' ;	'PV_' ;	'PW_' ;	'PY_' ;
'QA_' ;	'QC_' ;	'QD_' ;	'QE_' ;	'QF_' ;	'QG_' ;	'QH_' ;	'QI_' ;	'QK_' ;	'QL_' ;	'QM_' ;	'QN_' ;	'QP_' ;	'QQ_' ;	'QR_' ;	'QS_' ;	'QT_' ;	'QV_' ;	'QW_' ;	'QY_' ;
'RA_' ;	'RC_' ;	'RD_' ;	'RE_' ;	'RF_' ;	'RG_' ;	'RH_' ;	'RI_' ;	'RK_' ;	'RL_' ;	'RM_' ;	'RN_' ;	'RP_' ;	'RQ_' ;	'RR_' ;	'RS_' ;	'RT_' ;	'RV_' ;	'RW_' ;	'RY_' ;
'SA_' ;	'SC_' ;	'SD_' ;	'SE_' ;	'SF_' ;	'SG_' ;	'SH_' ;	'SI_' ;	'SK_' ;	'SL_' ;	'SM_' ;	'SN_' ;	'SP_' ;	'SQ_' ;	'SR_' ;	'SS_' ;	'ST_' ;	'SV_' ;	'SW_' ;	'SY_' ;
'TA_' ;	'TC_' ;	'TD_' ;	'TE_' ;	'TF_' ;	'TG_' ;	'TH_' ;	'TI_' ;	'TK_' ;	'TL_' ;	'TM_' ;	'TN_' ;	'TP_' ;	'TQ_' ;	'TR_' ;	'TS_' ;	'TT_' ;	'TV_' ;	'TW_' ;	'TY_' ;
'VA_' ;	'VC_' ;	'VD_' ;	'VE_' ;	'VF_' ;	'VG_' ;	'VH_' ;	'VI_' ;	'VK_' ;	'VL_' ;	'VM_' ;	'VN_' ;	'VP_' ;	'VQ_' ;	'VR_' ;	'VS_' ;	'VT_' ;	'VV_' ;	'VW_' ;	'VY_' ;
'WA_' ;	'WC_' ;	'WD_' ;	'WE_' ;	'WF_' ;	'WG_' ;	'WH_' ;	'WI_' ;	'WK_' ;	'WL_' ;	'WM_' ;	'WN_' ;	'WP_' ;	'WQ_' ;	'WR_' ;	'WS_' ;	'WT_' ;	'WV_' ;	'WW_' ;	'WY_' ;
'YA_' ;	'YC_' ;	'YD_' ;	'YE_' ;	'YF_' ;	'YG_' ;	'YH_' ;	'YI_' ;	'YK_' ;	'YL_' ;	'YM_' ;	'YN_' ;	'YP_' ;	'YQ_' ;	'YR_' ;	'YS_' ;	'YT_' ;	'YV_' ;	'YW_' ;	'YY_' ;
'A_A' ;	'A_C' ;	'A_D' ;	'A_E' ;	'A_F' ;	'A_G' ;	'A_H' ;	'A_I' ;	'A_K' ;	'A_L' ;	'A_M' ;	'A_N' ;	'A_P' ;	'A_Q' ;	'A_R' ;	'A_S' ;	'A_T' ;	'A_V' ;	'A_W' ;	'A_Y' ;
'C_A' ;	'C_C' ;	'C_D' ;	'C_E' ;	'C_F' ;	'C_G' ;	'C_H' ;	'C_I' ;	'C_K' ;	'C_L' ;	'C_M' ;	'C_N' ;	'C_P' ;	'C_Q' ;	'C_R' ;	'C_S' ;	'C_T' ;	'C_V' ;	'C_W' ;	'C_Y' ;
'D_A' ;	'D_C' ;	'D_D' ;	'D_E' ;	'D_F' ;	'D_G' ;	'D_H' ;	'D_I' ;	'D_K' ;	'D_L' ;	'D_M' ;	'D_N' ;	'D_P' ;	'D_Q' ;	'D_R' ;	'D_S' ;	'D_T' ;	'D_V' ;	'D_W' ;	'D_Y' ;
'E_A' ;	'E_C' ;	'E_D' ;	'E_E' ;	'E_F' ;	'E_G' ;	'E_H' ;	'E_I' ;	'E_K' ;	'E_L' ;	'E_M' ;	'E_N' ;	'E_P' ;	'E_Q' ;	'E_R' ;	'E_S' ;	'E_T' ;	'E_V' ;	'E_W' ;	'E_Y' ;
'F_A' ;	'F_C' ;	'F_D' ;	'F_E' ;	'F_F' ;	'F_G' ;	'F_H' ;	'F_I' ;	'F_K' ;	'F_L' ;	'F_M' ;	'F_N' ;	'F_P' ;	'F_Q' ;	'F_R' ;	'F_S' ;	'F_T' ;	'F_V' ;	'F_W' ;	'F_Y' ;
'G_A' ;	'G_C' ;	'G_D' ;	'G_E' ;	'G_F' ;	'G_G' ;	'G_H' ;	'G_I' ;	'G_K' ;	'G_L' ;	'G_M' ;	'G_N' ;	'G_P' ;	'G_Q' ;	'G_R' ;	'G_S' ;	'G_T' ;	'G_V' ;	'G_W' ;	'G_Y' ;
'H_A' ;	'H_C' ;	'H_D' ;	'H_E' ;	'H_F' ;	'H_G' ;	'H_H' ;	'H_I' ;	'H_K' ;	'H_L' ;	'H_M' ;	'H_N' ;	'H_P' ;	'H_Q' ;	'H_R' ;	'H_S' ;	'H_T' ;	'H_V' ;	'H_W' ;	'H_Y' ;
'I_A' ;	'I_C' ;	'I_D' ;	'I_E' ;	'I_F' ;	'I_G' ;	'I_H' ;	'I_I' ;	'I_K' ;	'I_L' ;	'I_M' ;	'I_N' ;	'I_P' ;	'I_Q' ;	'I_R' ;	'I_S' ;	'I_T' ;	'I_V' ;	'I_W' ;	'I_Y' ;
'K_A' ;	'K_C' ;	'K_D' ;	'K_E' ;	'K_F' ;	'K_G' ;	'K_H' ;	'K_I' ;	'K_K' ;	'K_L' ;	'K_M' ;	'K_N' ;	'K_P' ;	'K_Q' ;	'K_R' ;	'K_S' ;	'K_T' ;	'K_V' ;	'K_W' ;	'K_Y' ;
'L_A' ;	'L_C' ;	'L_D' ;	'L_E' ;	'L_F' ;	'L_G' ;	'L_H' ;	'L_I' ;	'L_K' ;	'L_L' ;	'L_M' ;	'L_N' ;	'L_P' ;	'L_Q' ;	'L_R' ;	'L_S' ;	'L_T' ;	'L_V' ;	'L_W' ;	'L_Y' ;
'M_A' ;	'M_C' ;	'M_D' ;	'M_E' ;	'M_F' ;	'M_G' ;	'M_H' ;	'M_I' ;	'M_K' ;	'M_L' ;	'M_M' ;	'M_N' ;	'M_P' ;	'M_Q' ;	'M_R' ;	'M_S' ;	'M_T' ;	'M_V' ;	'M_W' ;	'M_Y' ;
'N_A' ;	'N_C' ;	'N_D' ;	'N_E' ;	'N_F' ;	'N_G' ;	'N_H' ;	'N_I' ;	'N_K' ;	'N_L' ;	'N_M' ;	'N_N' ;	'N_P' ;	'N_Q' ;	'N_R' ;	'N_S' ;	'N_T' ;	'N_V' ;	'N_W' ;	'N_Y' ;
'P_A' ;	'P_C' ;	'P_D' ;	'P_E' ;	'P_F' ;	'P_G' ;	'P_H' ;	'P_I' ;	'P_K' ;	'P_L' ;	'P_M' ;	'P_N' ;	'P_P' ;	'P_Q' ;	'P_R' ;	'P_S' ;	'P_T' ;	'P_V' ;	'P_W' ;	'P_Y' ;
'Q_A' ;	'Q_C' ;	'Q_D' ;	'Q_E' ;	'Q_F' ;	'Q_G' ;	'Q_H' ;	'Q_I' ;	'Q_K' ;	'Q_L' ;	'Q_M' ;	'Q_N' ;	'Q_P' ;	'Q_Q' ;	'Q_R' ;	'Q_S' ;	'Q_T' ;	'Q_V' ;	'Q_W' ;	'Q_Y' ;
'R_A' ;	'R_C' ;	'R_D' ;	'R_E' ;	'R_F' ;	'R_G' ;	'R_H' ;	'R_I' ;	'R_K' ;	'R_L' ;	'R_M' ;	'R_N' ;	'R_P' ;	'R_Q' ;	'R_R' ;	'R_S' ;	'R_T' ;	'R_V' ;	'R_W' ;	'R_Y' ;
'S_A' ;	'S_C' ;	'S_D' ;	'S_E' ;	'S_F' ;	'S_G' ;	'S_H' ;	'S_I' ;	'S_K' ;	'S_L' ;	'S_M' ;	'S_N' ;	'S_P' ;	'S_Q' ;	'S_R' ;	'S_S' ;	'S_T' ;	'S_V' ;	'S_W' ;	'S_Y' ;
'T_A' ;	'T_C' ;	'T_D' ;	'T_E' ;	'T_F' ;	'T_G' ;	'T_H' ;	'T_I' ;	'T_K' ;	'T_L' ;	'T_M' ;	'T_N' ;	'T_P' ;	'T_Q' ;	'T_R' ;	'T_S' ;	'T_T' ;	'T_V' ;	'T_W' ;	'T_Y' ;
'V_A' ;	'V_C' ;	'V_D' ;	'V_E' ;	'V_F' ;	'V_G' ;	'V_H' ;	'V_I' ;	'V_K' ;	'V_L' ;	'V_M' ;	'V_N' ;	'V_P' ;	'V_Q' ;	'V_R' ;	'V_S' ;	'V_T' ;	'V_V' ;	'V_W' ;	'V_Y' ;
'W_A' ;	'W_C' ;	'W_D' ;	'W_E' ;	'W_F' ;	'W_G' ;	'W_H' ;	'W_I' ;	'W_K' ;	'W_L' ;	'W_M' ;	'W_N' ;	'W_P' ;	'W_Q' ;	'W_R' ;	'W_S' ;	'W_T' ;	'W_V' ;	'W_W' ;	'W_Y' ;
'Y_A' ;	'Y_C' ;	'Y_D' ;	'Y_E' ;	'Y_F' ;	'Y_G' ;	'Y_H' ;	'Y_I' ;	'Y_K' ;	'Y_L' ;	'Y_M' ;	'Y_N' ;	'Y_P' ;	'Y_Q' ;	'Y_R' ;	'Y_S' ;	'Y_T' ;	'Y_V' ;	'Y_W' ;	'Y_Y' ;
'X_A' ; 'X_C' ; 'X_D' ; 'X_E' ; 'X_F' ; 'X_G' ; 'X_H' ; 'X_I' ; 'X_K' ; 'X_L' ; 'X_M' ; 'X_N' ; 'X_P' ; 'X_Q' ; 'X_R' ; 'X_S' ; 'X_T' ; 'X_V' ; 'X_W' ; 'X_Y' ; 'XA_' ; 'XC_' ; 'XD_' ; 'XE_' ; 'XF_' ; 'XG_' ; 'XH_' ; 'XI_' ; 'XK_' ; 'XL_' ; 'XM_' ;
'XN_' ; 'XP_' ; 'XQ_' ; 'XR_' ; 'XS_' ; 'XT_' ; 'XV_' ; 'XX_' ; 'XW_' ; 'XY_' ; 'A_X' ; 'C_X' ; 'D_X' ; 'E_X' ; 'F_X' ; 'G_X' ; 'H_X' ; 'I_X' ; 'K_X' ; 'L_X' ; 'M_X' ; 'N_X' ; 'P_X' ; 'Q_X' ; 'R_X' ; 'S_X' ; 'T_X' ; 'V_X' ; 'X_X' ; 'W_X' ; 'Y_X' ; 'AX_' ;
'CX_' ; 'DX_' ; 'EX_' ; 'FX_' ; 'GX_' ; 'HX_' ; 'IX_' ; 'KX_' ; 'LX_' ; 'MX_' ; 'NX_' ; 'PX_' ; 'QX_' ; 'RX_' ; 'SX_' ; 'TX_' ; 'VX_' ; 'WX_' ; 'YX_' ;}; 

%% Simplified 3mer generation

Output_seq = zeros(length(TwomerPlus),length(Enrichment));
Output_count = zeros(1, length(Enrichment));

for n = 1: length(Enrichment);
    string_map1 = {};
    string_map2 = {};
    for a = 1: length(Enrichment{n})-1;
        string_map1{a,1} = cat(2, Enrichment{n}(a:a+1), '_');
    end
    for a = 1: length(Enrichment{n})-2;
        string_map2{a,1} = cat(2, Enrichment{n}(a), '_', Enrichment{n}(a+2));
    
    end
    string_map = cat(1, string_map1, string_map2);
    for x = 1: length(string_map);
        [kmer, row] = ismember(string_map(x), TwomerPlus);
        if kmer == 1;
            Output_seq(row, n) = Output_seq(row, n) + 1;
        else
            continue
        end
    end
end
% This routine breaks down each sequence on the Enrichment matrix into
% simplified 3mers and maps each 3mer to the TwomerPlus library. 
% In practice, it converts each sequence into a vector - in the case of a
% simplified 3mer library, this will be an 800-dimensional vector.


for n = 1 : size(Output_seq,2);
    Size_vec = 0;
    for x = 1 : size(Output_seq,1);
        Size_vec = Size_vec + Output_seq(x,n)^2;
    end
    Output_seq(:,n) = Output_seq(:,n)/sqrt(Size_vec);
end
% This routine converts all vectors to unit vectors

Output_count = cell2mat(Enrichment(:,2).');
% This aggregates the frequency of each sequence on a separate matrix

for n = 1 : length(Output_count);
    Output_seq(:,n) = Output_seq(:,n) * Output_count(n);
end
% This routine multiples the vector by its frequency, i.e. the vector
% length becomes a score of its enrichment. The impact of this on the PCA
% is still untested. The hypothesis is that this scaling introduces a bias
% in the PCA which should benefit more common sequences.


%% PCA

Output_seq = Output_seq.';
[Contrib_kmer,~,~,~,component,~] = pca(Output_seq);
% This generates two variables:
% 1. 'component' that lists the imapct of each dimension
% 2. 'Contrib_kmer' that provides the coefficiencts of the PCA for each
% 3mer

%% Picking out the more significant values v.3

cutoff = 0.05;


key_mer = {};
key_value = [];

for a = 1 : length(Contrib_kmer);          
  if Contrib_kmer(a,1) >= cutoff;
      if ismember(TwomerPlus(a), key_mer) == 1;
          key_value = key_value +  Contrib_kmer(a,1) * strcmp(key_mer, TwomerPlus(a));
      else
          key_mer = cat(1, key_mer, TwomerPlus(a));
          key_value = cat(1, key_value, 0);
          key_value = key_value +  Contrib_kmer(a,1) * strcmp(key_mer, TwomerPlus(a));
      end 
      
  elseif Contrib_kmer(a,1) <= -cutoff;
      if ismember(TwomerPlus(a), key_mer) == 1;
          key_value = key_value + Contrib_kmer(a,1) * strcmp(key_mer, TwomerPlus(a));
      else
          key_mer = cat(1, key_mer, TwomerPlus(a));
          key_value = cat(1, key_value, 0);
          key_value = key_value + Contrib_kmer(a,1) * strcmp(key_mer, TwomerPlus(a));
      end

  end
end

key_table1 = cat(2, key_mer, num2cell(key_value));   

key_mer = {};
key_value = [];

for a = 1 : length(Contrib_kmer);          
  if Contrib_kmer(a,2) >= cutoff;
      if ismember(TwomerPlus(a), key_mer) == 1;
          key_value = key_value +  Contrib_kmer(a,2) * strcmp(key_mer, TwomerPlus(a));
      else
          key_mer = cat(1, key_mer, TwomerPlus(a));
          key_value = cat(1, key_value, 0);
          key_value = key_value +  Contrib_kmer(a,2) * strcmp(key_mer, TwomerPlus(a));
      end 
      
  elseif Contrib_kmer(a,2) <= -cutoff;
      if ismember(TwomerPlus(a), key_mer) == 1;
          key_value = key_value + Contrib_kmer(a,2) * strcmp(key_mer, TwomerPlus(a));
      else
          key_mer = cat(1, key_mer, TwomerPlus(a));
          key_value = cat(1, key_value, 0);
          key_value = key_value + Contrib_kmer(a,2) * strcmp(key_mer, TwomerPlus(a));
      end

  end
end

key_table2 = cat(2, key_mer, num2cell(key_value));  

key_mer = {};
key_value = [];

for a = 1 : length(Contrib_kmer);          
  if Contrib_kmer(a,3) >= cutoff;
      if ismember(TwomerPlus(a), key_mer) == 1;
          key_value = key_value +  Contrib_kmer(a,3) * strcmp(key_mer, TwomerPlus(a));
      else
          key_mer = cat(1, key_mer, TwomerPlus(a));
          key_value = cat(1, key_value, 0);
          key_value = key_value +  Contrib_kmer(a,3) * strcmp(key_mer, TwomerPlus(a));
      end 
      
  elseif Contrib_kmer(a,3) <= - cutoff;
      if ismember(TwomerPlus(a), key_mer) == 1;
          key_value = key_value + Contrib_kmer(a,3) * strcmp(key_mer, TwomerPlus(a));
      else
          key_mer = cat(1, key_mer, TwomerPlus(a));
          key_value = cat(1, key_value, 0);
          key_value = key_value + Contrib_kmer(a,3) * strcmp(key_mer, TwomerPlus(a));
      end

  end
end

key_table3 = cat(2, key_mer, num2cell(key_value));  

key_mer = {};
key_value = [];

for a = 1 : length(Contrib_kmer);          
  if Contrib_kmer(a,4) >= cutoff;
      if ismember(TwomerPlus(a), key_mer) == 1;
          key_value = key_value +  Contrib_kmer(a,4) * strcmp(key_mer, TwomerPlus(a));
      else
          key_mer = cat(1, key_mer, TwomerPlus(a));
          key_value = cat(1, key_value, 0);
          key_value = key_value +  Contrib_kmer(a,4) * strcmp(key_mer, TwomerPlus(a));
      end 
      
  elseif Contrib_kmer(a,4) <= -cutoff;
      if ismember(TwomerPlus(a), key_mer) == 1;
          key_value = key_value + Contrib_kmer(a,4) * strcmp(key_mer, TwomerPlus(a));
      else
          key_mer = cat(1, key_mer, TwomerPlus(a));
          key_value = cat(1, key_value, 0);
          key_value = key_value + Contrib_kmer(a,4) * strcmp(key_mer, TwomerPlus(a));
      end

  end
end

key_table4 = cat(2, key_mer, num2cell(key_value));  

key_mer = {};
key_value = [];

for a = 1 : length(Contrib_kmer);          
  if Contrib_kmer(a,5) >= cutoff;
      if ismember(TwomerPlus(a), key_mer) == 1;
          key_value = key_value +  Contrib_kmer(a,5) * strcmp(key_mer, TwomerPlus(a));
      else
          key_mer = cat(1, key_mer, TwomerPlus(a));
          key_value = cat(1, key_value, 0);
          key_value = key_value +  Contrib_kmer(a,5) * strcmp(key_mer, TwomerPlus(a));
      end 
      
  elseif Contrib_kmer(a,5) <= -cutoff;
      if ismember(TwomerPlus(a), key_mer) == 1;
          key_value = key_value + Contrib_kmer(a,5) * strcmp(key_mer, TwomerPlus(a));
      else
          key_mer = cat(1, key_mer, TwomerPlus(a));
          key_value = cat(1, key_value, 0);
          key_value = key_value + Contrib_kmer(a,5) * strcmp(key_mer, TwomerPlus(a));
      end

  end
end

key_table5 = cat(2, key_mer, num2cell(key_value));  
                  

 contribution_CDF = [component(1), sum(component(1:2)), sum(component(1:3)), sum(component(1:4)), sum(component(1:5))];
 figure(5);
 bar(contribution_CDF, 1, 'FaceColor', [1,1,1], 'EdgeColor', [0,0,0], 'LineWidth', 1)
 
 
% This routine variation picks out all significant contributions whether
% they are positive or negative from the first 5 dimensions.

% varlist = {'key_mer', 'key_value', 'kmer', 'Output_count', 'Output_seq', 'Size_vec', 'string_map1', 'string_map2', 'x', 'varlist'};
% clear(varlist{:});