#!/usr/bin/env python
# encoding: utf-8

name = "H_Abstraction/TS_groups"
shortDesc = u""
longDesc = u"""

"""

entry(
    index = 0,
    label = "X_H_or_Xrad_H_Xbirad_H_Xtrirad_H",
    group = "OR{H2, C_H, O_H}",
    kinetics = DistanceData(
        distances = {'d12': 1.29714, 'd13': 2.57561, 'd23': 1.29365},
        uncertainties = {'d12': 0.162671, 'd13': 0.146329, 'd23': 0.146651},
    ),
    shortDesc = u"""Fitted to 2490 distances.
""",
    longDesc = 
u"""
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=166 label="OjH">]
[<Entry index=2 label="C_H">, <Entry index=168 label="OjCs">]
[<Entry index=155 label="OHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=151 label="Cb_H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=151 label="Cb_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=268 label="CtjC">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=150 label="Ct_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=155 label="OHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=2 label="C_H">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=269 label="Cbj">]
[<Entry index=5 label="C_methane">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=1 label="H2">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=4 label="Csnorad_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=2 label="C_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=153 label="OradH">, <Entry index=163 label="Hrad">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=155 label="OHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=166 label="OjH">]
[<Entry index=5 label="C_methane">, <Entry index=166 label="OjH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=269 label="Cbj">]
[<Entry index=7 label="CsCHHH">, <Entry index=166 label="OjH">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=25 label="CsCOHH">, <Entry index=172 label="OjO">]
[<Entry index=151 label="Cb_H">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=2 label="C_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=155 label="OHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=5 label="C_methane">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=7 label="CsCHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=2 label="C_H">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=155 label="OHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=71 label="C_methyl">, <Entry index=163 label="Hrad">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=5 label="C_methane">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=166 label="OjH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=161 label="OOH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=2 label="C_H">, <Entry index=172 label="OjO">]
[<Entry index=155 label="OHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=1 label="H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=163 label="Hrad">]
[<Entry index=7 label="CsCHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=71 label="C_methyl">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=5 label="C_methane">, <Entry index=175 label="Cj">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=161 label="OOH">, <Entry index=186 label="CsjCCH">]
[<Entry index=151 label="Cb_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=163 label="Hrad">]
[<Entry index=71 label="C_methyl">, <Entry index=204 label="CsjCCC">]
[<Entry index=71 label="C_methyl">, <Entry index=172 label="OjO">]
[<Entry index=12 label="CsOHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=204 label="CsjCCC">]
[<Entry index=12 label="CsOHHH">, <Entry index=269 label="Cbj">]
[<Entry index=161 label="OOH">, <Entry index=163 label="Hrad">]
[<Entry index=71 label="C_methyl">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=1 label="H2">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=153 label="OradH">, <Entry index=186 label="CsjCCH">]
[<Entry index=1 label="H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=161 label="OOH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=5 label="C_methane">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=12 label="CsOHHH">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=161 label="OOH">, <Entry index=184 label="CsjOH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=175 label="Cj">]
[<Entry index=161 label="OOH">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=71 label="C_methyl">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=7 label="CsCHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=7 label="CsCHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=153 label="OradH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=1 label="H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=175 label="Cj">]
[<Entry index=1 label="H2">, <Entry index=172 label="OjO">]
[<Entry index=154 label="ORH">, <Entry index=184 label="CsjOH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=184 label="CsjOH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=153 label="OradH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=151 label="Cb_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=177 label="Cs_methyl">]
[<Entry index=155 label="OHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=1 label="H2">, <Entry index=166 label="OjH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=150 label="Ct_H">, <Entry index=268 label="CtjC">]
[<Entry index=150 label="Ct_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=175 label="Cj">]
[<Entry index=150 label="Ct_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=1 label="H2">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=172 label="OjO">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=1 label="H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=150 label="Ct_H">, <Entry index=269 label="Cbj">]
[<Entry index=4 label="Csnorad_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=1 label="H2">, <Entry index=175 label="Cj">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=150 label="Ct_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=151 label="Cb_H">, <Entry index=163 label="Hrad">]
[<Entry index=154 label="ORH">, <Entry index=172 label="OjO">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=163 label="Hrad">]
[<Entry index=5 label="C_methane">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=179 label="CsjCH2">]
[<Entry index=127 label="Cd_Cds/Cd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=2 label="C_H">, <Entry index=175 label="Cj">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=204 label="CsjCCC">]
[<Entry index=150 label="Ct_H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=71 label="C_methyl">, <Entry index=175 label="Cj">]
[<Entry index=161 label="OOH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=174 label="Crad">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=161 label="OOH">, <Entry index=179 label="CsjCH2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=166 label="OjH">]
[<Entry index=161 label="OOH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=150 label="Ct_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=153 label="OradH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=1 label="H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=172 label="OjO">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=4 label="Csnorad_H">, <Entry index=168 label="OjCs">]
[<Entry index=7 label="CsCHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=161 label="OOH">, <Entry index=175 label="Cj">]
[<Entry index=161 label="OOH">, <Entry index=172 label="OjO">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=150 label="Ct_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=155 label="OHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=155 label="OHH">, <Entry index=172 label="OjO">]
[<Entry index=153 label="OradH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=5 label="C_methane">, <Entry index=172 label="OjO">]
[<Entry index=151 label="Cb_H">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=174 label="Crad">]
[<Entry index=12 label="CsOHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=175 label="Cj">]
[<Entry index=14 label="CsCCHH">, <Entry index=166 label="OjH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=175 label="Cj">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=175 label="Cj">]
[<Entry index=155 label="OHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=1 label="H2">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=12 label="CsOHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=1 label="H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=168 label="OjCs">]
[<Entry index=153 label="OradH">, <Entry index=204 label="CsjCCC">]
[<Entry index=14 label="CsCCHH">, <Entry index=172 label="OjO">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=172 label="OjO">]
[<Entry index=153 label="OradH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=268 label="CtjC">]
[<Entry index=121 label="Cd_C/R/H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=163 label="Hrad">]
[<Entry index=12 label="CsOHHH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=151 label="Cb_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=153 label="OradH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=1 label="H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=5 label="C_methane">, <Entry index=186 label="CsjCCH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=4 label="Csnorad_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=163 label="Hrad">]
[<Entry index=151 label="Cb_H">, <Entry index=175 label="Cj">]
[<Entry index=161 label="OOH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=151 label="Cb_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=7 label="CsCHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=175 label="Cj">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=174 label="Crad">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=172 label="OjO">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=155 label="OHH">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=161 label="OOH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=2 label="C_H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=2 label="C_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=7 label="CsCHHH">, <Entry index=168 label="OjCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=71 label="C_methyl">, <Entry index=166 label="OjH">]
[<Entry index=7 label="CsCHHH">, <Entry index=172 label="OjO">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=172 label="OjO">]
[<Entry index=2 label="C_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=153 label="OradH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=2 label="C_H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=166 label="OjH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=151 label="Cb_H">, <Entry index=172 label="OjO">]
[<Entry index=1 label="H2">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=5 label="C_methane">, <Entry index=168 label="OjCs">]
[<Entry index=153 label="OradH">, <Entry index=175 label="Cj">]
[<Entry index=153 label="OradH">, <Entry index=179 label="CsjCH2">]
[<Entry index=5 label="C_methane">, <Entry index=169 label="OjCd">]
[<Entry index=151 label="Cb_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=153 label="OradH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=12 label="CsOHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=4 label="Csnorad_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=151 label="Cb_H">, <Entry index=269 label="Cbj">]
[<Entry index=155 label="OHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=204 label="CsjCCC">]
[<Entry index=5 label="C_methane">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=166 label="OjH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=186 label="CsjCCH">]
[<Entry index=1 label="H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=168 label="OjCs">]
[<Entry index=1 label="H2">, <Entry index=169 label="OjCd">]
[<Entry index=71 label="C_methyl">, <Entry index=168 label="OjCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=5 label="C_methane">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=153 label="OradH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=1 label="H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=5 label="C_methane">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=151 label="Cb_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=161 label="OOH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=151 label="Cb_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=172 label="OjO">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=204 label="CsjCCC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=168 label="OjCs">]
[<Entry index=151 label="Cb_H">, <Entry index=166 label="OjH">]
[<Entry index=1 label="H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=166 label="OjH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=172 label="OjO">]
[<Entry index=4 label="Csnorad_H">, <Entry index=175 label="Cj">]
[<Entry index=12 label="CsOHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=163 label="Hrad">]
[<Entry index=161 label="OOH">, <Entry index=169 label="OjCd">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=179 label="CsjCH2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=184 label="CsjOH2">]
[<Entry index=153 label="OradH">, <Entry index=168 label="OjCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=268 label="CtjC">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=1 label="H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=161 label="OOH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=2 label="C_H">, <Entry index=166 label="OjH">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=172 label="OjO">]
[<Entry index=1 label="H2">, <Entry index=204 label="CsjCCC">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=5 label="C_methane">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=12 label="CsOHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=186 label="CsjCCH">]
[<Entry index=5 label="C_methane">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=14 label="CsCCHH">, <Entry index=163 label="Hrad">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=155 label="OHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=155 label="OHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=2 label="C_H">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=269 label="Cbj">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=172 label="OjO">]
[<Entry index=161 label="OOH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=168 label="OjCs">]
[<Entry index=2 label="C_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=1 label="H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=7 label="CsCHHH">, <Entry index=163 label="Hrad">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=5 label="C_methane">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=172 label="OjO">]
[<Entry index=155 label="OHH">, <Entry index=335 label="C_quartetR">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=161 label="OOH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=14 label="CsCCHH">, <Entry index=168 label="OjCs">]
[<Entry index=71 label="C_methyl">, <Entry index=186 label="CsjCCH">]
[<Entry index=25 label="CsCOHH">, <Entry index=163 label="Hrad">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=7 label="CsCHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=1 label="H2">, <Entry index=168 label="OjCs">]
[<Entry index=1 label="H2">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=5 label="C_methane">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=71 label="C_methyl">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=155 label="OHH">, <Entry index=269 label="Cbj">]
[<Entry index=155 label="OHH">, <Entry index=168 label="OjCs">]
[<Entry index=150 label="Ct_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=7 label="CsCHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=5 label="C_methane">, <Entry index=163 label="Hrad">]
[<Entry index=1 label="H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=2 label="C_H">, <Entry index=269 label="Cbj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=268 label="CtjC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=163 label="Hrad">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=161 label="OOH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=2 label="C_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=14 label="CsCCHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=161 label="OOH">, <Entry index=168 label="OjCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=269 label="Cbj">]
[<Entry index=150 label="Ct_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=161 label="OOH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=175 label="Cj">]
[<Entry index=4 label="Csnorad_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=153 label="OradH">, <Entry index=184 label="CsjOH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=2 label="C_H">, <Entry index=167 label="OjC">]
[<Entry index=155 label="OHH">, <Entry index=175 label="Cj">]
[<Entry index=153 label="OradH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=151 label="Cb_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=1 label="H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=4 label="Csnorad_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=12 label="CsOHHH">, <Entry index=166 label="OjH">]
[<Entry index=155 label="OHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=151 label="Cb_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=172 label="OjO">]
[<Entry index=2 label="C_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=153 label="OradH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=150 label="Ct_H">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=269 label="Cbj">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=153 label="OradH">, <Entry index=172 label="OjO">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=1 label="H2">, <Entry index=186 label="CsjCCH">]
[<Entry index=1 label="H2">, <Entry index=335 label="C_quartetR">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=5 label="C_methane">, <Entry index=184 label="CsjOH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=163 label="Hrad">]
[<Entry index=153 label="OradH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=166 label="OjH">]
[<Entry index=155 label="OHH">, <Entry index=163 label="Hrad">]
[<Entry index=2 label="C_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=5 label="C_methane">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=5 label="C_methane">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=12 label="CsOHHH">, <Entry index=167 label="OjC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=186 label="CsjCCH">]
[<Entry index=7 label="CsCHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=163 label="Hrad">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=15 label="C/H2/Cs/Cs">, <Entry index=172 label="OjO">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=177 label="Cs_methyl">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=7 label="CsCHHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=150 label="Ct_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=5 label="C_methane">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=2 label="C_H">, <Entry index=163 label="Hrad">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=168 label="OjCs">]
[<Entry index=5 label="C_methane">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=2 label="C_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=2 label="C_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=166 label="OjH">]
[<Entry index=151 label="Cb_H">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=1 label="H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=25 label="CsCOHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=12 label="CsOHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=2 label="C_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=166 label="OjH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=335 label="C_quartetR">]
[<Entry index=7 label="CsCHHH">, <Entry index=269 label="Cbj">]
""",
)

entry(
    index = 1,
    label = "Y_rad_birad_trirad_quadrad",
    group = "OR{Hrad, Orad, Crad}",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 2,
    label = "H2",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 H u0 {2,S}
2 *2 H u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.332248, 'd13': -0.344377, 'd23': -0.023575},
        uncertainties = {'d12': 0.081855, 'd13': 0.047408, 'd23': 0.071279},
    ),
    shortDesc = u"""Fitted to 101 distances.
""",
    longDesc = 
u"""
[<Entry index=1 label="H2">, <Entry index=168 label="OjCs">]
[<Entry index=1 label="H2">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=1 label="H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=1 label="H2">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=1 label="H2">, <Entry index=175 label="Cj">]
[<Entry index=1 label="H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=1 label="H2">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=1 label="H2">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=1 label="H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=1 label="H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=1 label="H2">, <Entry index=186 label="CsjCCH">]
[<Entry index=1 label="H2">, <Entry index=335 label="C_quartetR">]
[<Entry index=1 label="H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=1 label="H2">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=1 label="H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=1 label="H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=1 label="H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=1 label="H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=1 label="H2">, <Entry index=172 label="OjO">]
[<Entry index=1 label="H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=1 label="H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=1 label="H2">, <Entry index=169 label="OjCd">]
[<Entry index=1 label="H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=1 label="H2">, <Entry index=166 label="OjH">]
[<Entry index=1 label="H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=1 label="H2">, <Entry index=204 label="CsjCCC">]
[<Entry index=1 label="H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=1 label="H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=1 label="H2">, <Entry index=173 label="O_atom_triplet">]
""",
)

entry(
    index = 3,
    label = "C_H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 C ux {2,S}
2 *2 H u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.021334, 'd13': 0.028247, 'd23': 0.00654},
        uncertainties = {'d12': 0.167859, 'd13': 0.149632, 'd23': 0.146199},
    ),
    shortDesc = u"""Fitted to 2145 distances.
""",
    longDesc = 
u"""
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=166 label="OjH">]
[<Entry index=2 label="C_H">, <Entry index=168 label="OjCs">]
[<Entry index=151 label="Cb_H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=151 label="Cb_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=268 label="CtjC">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=150 label="Ct_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=2 label="C_H">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=269 label="Cbj">]
[<Entry index=5 label="C_methane">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=184 label="CsjOH2">]
[<Entry index=4 label="Csnorad_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=2 label="C_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=166 label="OjH">]
[<Entry index=5 label="C_methane">, <Entry index=166 label="OjH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=269 label="Cbj">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=204 label="CsjCCC">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=25 label="CsCOHH">, <Entry index=172 label="OjO">]
[<Entry index=151 label="Cb_H">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=2 label="C_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=5 label="C_methane">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=7 label="CsCHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=2 label="C_H">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=71 label="C_methyl">, <Entry index=163 label="Hrad">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=5 label="C_methane">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=166 label="OjH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=2 label="C_H">, <Entry index=172 label="OjO">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=5 label="C_methane">, <Entry index=163 label="Hrad">]
[<Entry index=7 label="CsCHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=71 label="C_methyl">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=5 label="C_methane">, <Entry index=175 label="Cj">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=151 label="Cb_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=163 label="Hrad">]
[<Entry index=71 label="C_methyl">, <Entry index=204 label="CsjCCC">]
[<Entry index=71 label="C_methyl">, <Entry index=172 label="OjO">]
[<Entry index=12 label="CsOHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=204 label="CsjCCC">]
[<Entry index=12 label="CsOHHH">, <Entry index=269 label="Cbj">]
[<Entry index=71 label="C_methyl">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=5 label="C_methane">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=12 label="CsOHHH">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=7 label="CsCHHH">, <Entry index=175 label="Cj">]
[<Entry index=71 label="C_methyl">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=7 label="CsCHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=7 label="CsCHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=175 label="Cj">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=151 label="Cb_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=177 label="Cs_methyl">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=150 label="Ct_H">, <Entry index=268 label="CtjC">]
[<Entry index=150 label="Ct_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=175 label="Cj">]
[<Entry index=150 label="Ct_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=172 label="OjO">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=150 label="Ct_H">, <Entry index=269 label="Cbj">]
[<Entry index=7 label="CsCHHH">, <Entry index=168 label="OjCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=174 label="Crad">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=150 label="Ct_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=151 label="Cb_H">, <Entry index=163 label="Hrad">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=163 label="Hrad">]
[<Entry index=5 label="C_methane">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=127 label="Cd_Cds/Cd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=2 label="C_H">, <Entry index=175 label="Cj">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=204 label="CsjCCC">]
[<Entry index=150 label="Ct_H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=174 label="Crad">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=166 label="OjH">]
[<Entry index=150 label="Ct_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=172 label="OjO">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=4 label="Csnorad_H">, <Entry index=168 label="OjCs">]
[<Entry index=7 label="CsCHHH">, <Entry index=166 label="OjH">]
[<Entry index=7 label="CsCHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=71 label="C_methyl">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=5 label="C_methane">, <Entry index=172 label="OjO">]
[<Entry index=151 label="Cb_H">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=175 label="Cj">]
[<Entry index=14 label="CsCCHH">, <Entry index=166 label="OjH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=175 label="Cj">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=175 label="Cj">]
[<Entry index=12 label="CsOHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=168 label="OjCs">]
[<Entry index=151 label="Cb_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=172 label="OjO">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=268 label="CtjC">]
[<Entry index=121 label="Cd_C/R/H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=163 label="Hrad">]
[<Entry index=12 label="CsOHHH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=151 label="Cb_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=5 label="C_methane">, <Entry index=186 label="CsjCCH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=4 label="Csnorad_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=163 label="Hrad">]
[<Entry index=151 label="Cb_H">, <Entry index=175 label="Cj">]
[<Entry index=151 label="Cb_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=7 label="CsCHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=175 label="Cj">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=174 label="Crad">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=172 label="OjO">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=2 label="C_H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=2 label="C_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=4 label="Csnorad_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=71 label="C_methyl">, <Entry index=166 label="OjH">]
[<Entry index=7 label="CsCHHH">, <Entry index=172 label="OjO">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=172 label="OjO">]
[<Entry index=2 label="C_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=2 label="C_H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=166 label="OjH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=151 label="Cb_H">, <Entry index=172 label="OjO">]
[<Entry index=5 label="C_methane">, <Entry index=168 label="OjCs">]
[<Entry index=5 label="C_methane">, <Entry index=169 label="OjCd">]
[<Entry index=151 label="Cb_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=12 label="CsOHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=4 label="Csnorad_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=151 label="Cb_H">, <Entry index=269 label="Cbj">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=204 label="CsjCCC">]
[<Entry index=5 label="C_methane">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=166 label="OjH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=186 label="CsjCCH">]
[<Entry index=12 label="CsOHHH">, <Entry index=168 label="OjCs">]
[<Entry index=71 label="C_methyl">, <Entry index=168 label="OjCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=5 label="C_methane">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=5 label="C_methane">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=151 label="Cb_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=151 label="Cb_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=172 label="OjO">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=168 label="OjCs">]
[<Entry index=151 label="Cb_H">, <Entry index=166 label="OjH">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=166 label="OjH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=172 label="OjO">]
[<Entry index=4 label="Csnorad_H">, <Entry index=175 label="Cj">]
[<Entry index=12 label="CsOHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=163 label="Hrad">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=179 label="CsjCH2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=184 label="CsjOH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=268 label="CtjC">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=2 label="C_H">, <Entry index=166 label="OjH">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=172 label="OjO">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=5 label="C_methane">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=12 label="CsOHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=186 label="CsjCCH">]
[<Entry index=5 label="C_methane">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=14 label="CsCCHH">, <Entry index=163 label="Hrad">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=268 label="CtjC">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=2 label="C_H">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=269 label="Cbj">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=172 label="OjO">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=168 label="OjCs">]
[<Entry index=2 label="C_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=5 label="C_methane">, <Entry index=179 label="CsjCH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=163 label="Hrad">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=5 label="C_methane">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=172 label="OjO">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=14 label="CsCCHH">, <Entry index=168 label="OjCs">]
[<Entry index=71 label="C_methyl">, <Entry index=186 label="CsjCCH">]
[<Entry index=25 label="CsCOHH">, <Entry index=163 label="Hrad">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=5 label="C_methane">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=71 label="C_methyl">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=7 label="CsCHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=2 label="C_H">, <Entry index=269 label="Cbj">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=163 label="Hrad">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=2 label="C_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=14 label="CsCCHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=269 label="Cbj">]
[<Entry index=150 label="Ct_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=175 label="Cj">]
[<Entry index=4 label="Csnorad_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=2 label="C_H">, <Entry index=167 label="OjC">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=151 label="Cb_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=4 label="Csnorad_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=12 label="CsOHHH">, <Entry index=166 label="OjH">]
[<Entry index=14 label="CsCCHH">, <Entry index=172 label="OjO">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=172 label="OjO">]
[<Entry index=2 label="C_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=150 label="Ct_H">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=269 label="Cbj">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=5 label="C_methane">, <Entry index=184 label="CsjOH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=163 label="Hrad">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=166 label="OjH">]
[<Entry index=2 label="C_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=5 label="C_methane">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=163 label="Hrad">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=5 label="C_methane">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=12 label="CsOHHH">, <Entry index=167 label="OjC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=186 label="CsjCCH">]
[<Entry index=7 label="CsCHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=163 label="Hrad">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=15 label="C/H2/Cs/Cs">, <Entry index=172 label="OjO">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=177 label="Cs_methyl">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=7 label="CsCHHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=150 label="Ct_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=5 label="C_methane">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=2 label="C_H">, <Entry index=163 label="Hrad">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=168 label="OjCs">]
[<Entry index=5 label="C_methane">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=2 label="C_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=2 label="C_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=166 label="OjH">]
[<Entry index=151 label="Cb_H">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=25 label="CsCOHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=12 label="CsOHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=2 label="C_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=166 label="OjH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=335 label="C_quartetR">]
[<Entry index=7 label="CsCHHH">, <Entry index=269 label="Cbj">]
""",
)

entry(
    index = 4,
    label = "Cs_H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs ux {2,S}
2 *2 H  u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.006936, 'd13': 0.019039, 'd23': 0.009985},
        uncertainties = {'d12': 0.158766, 'd13': 0.152084, 'd23': 0.139634},
    ),
    shortDesc = u"""Fitted to 1719 distances.
""",
    longDesc = 
u"""
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=12 label="CsOHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=4 label="Csnorad_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=5 label="C_methane">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=71 label="C_methyl">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=166 label="OjH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=175 label="Cj">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=172 label="OjO">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=172 label="OjO">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=5 label="C_methane">, <Entry index=186 label="CsjCCH">]
[<Entry index=7 label="CsCHHH">, <Entry index=268 label="CtjC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=4 label="Csnorad_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=177 label="Cs_methyl">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=7 label="CsCHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=268 label="CtjC">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=5 label="C_methane">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=163 label="Hrad">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=7 label="CsCHHH">, <Entry index=168 label="OjCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=71 label="C_methyl">, <Entry index=166 label="OjH">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=172 label="OjO">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=168 label="OjCs">]
[<Entry index=71 label="C_methyl">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=14 label="CsCCHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=163 label="Hrad">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=71 label="C_methyl">, <Entry index=175 label="Cj">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=25 label="CsCOHH">, <Entry index=172 label="OjO">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=166 label="OjH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=5 label="C_methane">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=5 label="C_methane">, <Entry index=166 label="OjH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=269 label="Cbj">]
[<Entry index=7 label="CsCHHH">, <Entry index=166 label="OjH">]
[<Entry index=71 label="C_methyl">, <Entry index=172 label="OjO">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=4 label="Csnorad_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=5 label="C_methane">, <Entry index=168 label="OjCs">]
[<Entry index=5 label="C_methane">, <Entry index=179 label="CsjCH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=5 label="C_methane">, <Entry index=169 label="OjCd">]
[<Entry index=7 label="CsCHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=168 label="OjCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=12 label="CsOHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=5 label="C_methane">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=4 label="Csnorad_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=71 label="C_methyl">, <Entry index=186 label="CsjCCH">]
[<Entry index=12 label="CsOHHH">, <Entry index=166 label="OjH">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=5 label="C_methane">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=4 label="Csnorad_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=7 label="CsCHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=186 label="CsjCCH">]
[<Entry index=12 label="CsOHHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=5 label="C_methane">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=71 label="C_methyl">, <Entry index=163 label="Hrad">]
[<Entry index=12 label="CsOHHH">, <Entry index=168 label="OjCs">]
[<Entry index=71 label="C_methyl">, <Entry index=168 label="OjCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=5 label="C_methane">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=5 label="C_methane">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=174 label="Crad">]
[<Entry index=12 label="CsOHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=166 label="OjH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=166 label="OjH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=269 label="Cbj">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=163 label="Hrad">]
[<Entry index=5 label="C_methane">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=166 label="OjH">]
[<Entry index=7 label="CsCHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=5 label="C_methane">, <Entry index=163 label="Hrad">]
[<Entry index=7 label="CsCHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=5 label="C_methane">, <Entry index=175 label="Cj">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=172 label="OjO">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=12 label="CsOHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=172 label="OjO">]
[<Entry index=71 label="C_methyl">, <Entry index=204 label="CsjCCC">]
[<Entry index=4 label="Csnorad_H">, <Entry index=175 label="Cj">]
[<Entry index=12 label="CsOHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=184 label="CsjOH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=177 label="Cs_methyl">]
[<Entry index=4 label="Csnorad_H">, <Entry index=168 label="OjCs">]
[<Entry index=5 label="C_methane">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=163 label="Hrad">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=204 label="CsjCCC">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=166 label="OjH">]
[<Entry index=12 label="CsOHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=179 label="CsjCH2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=184 label="CsjOH2">]
[<Entry index=5 label="C_methane">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=268 label="CtjC">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=204 label="CsjCCC">]
[<Entry index=12 label="CsOHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=12 label="CsOHHH">, <Entry index=269 label="Cbj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=71 label="C_methyl">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=163 label="Hrad">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=7 label="CsCHHH">, <Entry index=172 label="OjO">]
[<Entry index=5 label="C_methane">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=172 label="OjO">]
[<Entry index=5 label="C_methane">, <Entry index=172 label="OjO">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=174 label="Crad">]
[<Entry index=12 label="CsOHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=167 label="OjC">]
[<Entry index=12 label="CsOHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=175 label="Cj">]
[<Entry index=7 label="CsCHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=5 label="C_methane">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=12 label="CsOHHH">, <Entry index=175 label="Cj">]
[<Entry index=5 label="C_methane">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=14 label="CsCCHH">, <Entry index=163 label="Hrad">]
[<Entry index=12 label="CsOHHH">, <Entry index=163 label="Hrad">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=14 label="CsCCHH">, <Entry index=166 label="OjH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=175 label="Cj">]
[<Entry index=15 label="C/H2/Cs/Cs">, <Entry index=172 label="OjO">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=175 label="Cj">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=335 label="C_quartetR">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=71 label="C_methyl">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=12 label="CsOHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=269 label="Cbj">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=5 label="C_methane">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=5 label="C_methane">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=7 label="CsCHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=7 label="CsCHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=163 label="Hrad">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=5 label="C_methane">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=14 label="CsCCHH">, <Entry index=172 label="OjO">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=172 label="OjO">]
[<Entry index=12 label="CsOHHH">, <Entry index=172 label="OjO">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=166 label="OjH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=4 label="Csnorad_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=163 label="Hrad">]
[<Entry index=14 label="CsCCHH">, <Entry index=168 label="OjCs">]
[<Entry index=12 label="CsOHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=7 label="CsCHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=12 label="CsOHHH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=25 label="CsCOHH">, <Entry index=163 label="Hrad">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=7 label="CsCHHH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=166 label="OjH">]
[<Entry index=25 label="CsCOHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=269 label="Cbj">]
[<Entry index=12 label="CsOHHH">, <Entry index=181 label="Csj/Cd/H2">]
""",
)

entry(
    index = 5,
    label = "Csnorad_H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S}
2 *2 H  u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.005815, 'd13': 0.018819, 'd23': 0.010924},
        uncertainties = {'d12': 0.159562, 'd13': 0.153033, 'd23': 0.140653},
    ),
    shortDesc = u"""Fitted to 1689 distances.
""",
    longDesc = 
u"""
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=12 label="CsOHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=4 label="Csnorad_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=5 label="C_methane">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=7 label="CsCHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=166 label="OjH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=175 label="Cj">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=172 label="OjO">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=172 label="OjO">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=5 label="C_methane">, <Entry index=186 label="CsjCCH">]
[<Entry index=7 label="CsCHHH">, <Entry index=268 label="CtjC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=4 label="Csnorad_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=177 label="Cs_methyl">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=268 label="CtjC">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=5 label="C_methane">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=163 label="Hrad">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=7 label="CsCHHH">, <Entry index=168 label="OjCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=174 label="Crad">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=172 label="OjO">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=168 label="OjCs">]
[<Entry index=14 label="CsCCHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=163 label="Hrad">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=25 label="CsCOHH">, <Entry index=172 label="OjO">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=166 label="OjH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=5 label="C_methane">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=5 label="C_methane">, <Entry index=166 label="OjH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=269 label="Cbj">]
[<Entry index=7 label="CsCHHH">, <Entry index=166 label="OjH">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=4 label="Csnorad_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=5 label="C_methane">, <Entry index=168 label="OjCs">]
[<Entry index=5 label="C_methane">, <Entry index=179 label="CsjCH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=5 label="C_methane">, <Entry index=169 label="OjCd">]
[<Entry index=7 label="CsCHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=168 label="OjCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=12 label="CsOHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=5 label="C_methane">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=4 label="Csnorad_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=12 label="CsOHHH">, <Entry index=166 label="OjH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=5 label="C_methane">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=4 label="Csnorad_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=7 label="CsCHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=186 label="CsjCCH">]
[<Entry index=12 label="CsOHHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=5 label="C_methane">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=5 label="C_methane">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=12 label="CsOHHH">, <Entry index=168 label="OjCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=5 label="C_methane">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=174 label="Crad">]
[<Entry index=12 label="CsOHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=166 label="OjH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=5 label="C_methane">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=269 label="Cbj">]
[<Entry index=5 label="C_methane">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=166 label="OjH">]
[<Entry index=7 label="CsCHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=5 label="C_methane">, <Entry index=163 label="Hrad">]
[<Entry index=7 label="CsCHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=5 label="C_methane">, <Entry index=175 label="Cj">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=172 label="OjO">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=12 label="CsOHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=172 label="OjO">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=4 label="Csnorad_H">, <Entry index=175 label="Cj">]
[<Entry index=12 label="CsOHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=184 label="CsjOH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=177 label="Cs_methyl">]
[<Entry index=4 label="Csnorad_H">, <Entry index=168 label="OjCs">]
[<Entry index=5 label="C_methane">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=163 label="Hrad">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=204 label="CsjCCC">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=166 label="OjH">]
[<Entry index=12 label="CsOHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=179 label="CsjCH2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=184 label="CsjOH2">]
[<Entry index=5 label="C_methane">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=268 label="CtjC">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=204 label="CsjCCC">]
[<Entry index=12 label="CsOHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=12 label="CsOHHH">, <Entry index=269 label="Cbj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=163 label="Hrad">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=7 label="CsCHHH">, <Entry index=172 label="OjO">]
[<Entry index=5 label="C_methane">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=172 label="OjO">]
[<Entry index=5 label="C_methane">, <Entry index=172 label="OjO">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=167 label="OjC">]
[<Entry index=12 label="CsOHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=5 label="C_methane">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=12 label="CsOHHH">, <Entry index=175 label="Cj">]
[<Entry index=5 label="C_methane">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=14 label="CsCCHH">, <Entry index=163 label="Hrad">]
[<Entry index=12 label="CsOHHH">, <Entry index=163 label="Hrad">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=14 label="CsCCHH">, <Entry index=166 label="OjH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=175 label="Cj">]
[<Entry index=15 label="C/H2/Cs/Cs">, <Entry index=172 label="OjO">]
[<Entry index=7 label="CsCHHH">, <Entry index=175 label="Cj">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=175 label="Cj">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=12 label="CsOHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=5 label="C_methane">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=7 label="CsCHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=7 label="CsCHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=163 label="Hrad">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=5 label="C_methane">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=14 label="CsCCHH">, <Entry index=172 label="OjO">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=172 label="OjO">]
[<Entry index=12 label="CsOHHH">, <Entry index=172 label="OjO">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=166 label="OjH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=4 label="Csnorad_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=335 label="C_quartetR">]
[<Entry index=14 label="CsCCHH">, <Entry index=168 label="OjCs">]
[<Entry index=12 label="CsOHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=7 label="CsCHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=12 label="CsOHHH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=25 label="CsCOHH">, <Entry index=163 label="Hrad">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=7 label="CsCHHH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=25 label="CsCOHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=269 label="Cbj">]
[<Entry index=12 label="CsOHHH">, <Entry index=181 label="Csj/Cd/H2">]
""",
)

entry(
    index = 6,
    label = "C_methane",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
5    H  u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.096022, 'd13': 0.064801, 'd23': -0.04245},
        uncertainties = {'d12': 0.082851, 'd13': 0.060239, 'd23': 0.04707},
    ),
    shortDesc = u"""Fitted to 78 distances.
""",
    longDesc = 
u"""
[<Entry index=5 label="C_methane">, <Entry index=186 label="CsjCCH">]
[<Entry index=5 label="C_methane">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=5 label="C_methane">, <Entry index=163 label="Hrad">]
[<Entry index=5 label="C_methane">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=5 label="C_methane">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=5 label="C_methane">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=5 label="C_methane">, <Entry index=168 label="OjCs">]
[<Entry index=5 label="C_methane">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=5 label="C_methane">, <Entry index=175 label="Cj">]
[<Entry index=5 label="C_methane">, <Entry index=169 label="OjCd">]
[<Entry index=5 label="C_methane">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=5 label="C_methane">, <Entry index=184 label="CsjOH2">]
[<Entry index=5 label="C_methane">, <Entry index=166 label="OjH">]
[<Entry index=5 label="C_methane">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=5 label="C_methane">, <Entry index=172 label="OjO">]
[<Entry index=5 label="C_methane">, <Entry index=179 label="CsjCH2">]
[<Entry index=5 label="C_methane">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=5 label="C_methane">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=5 label="C_methane">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=5 label="C_methane">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=5 label="C_methane">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=5 label="C_methane">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=5 label="C_methane">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=5 label="C_methane">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=250 label="Cdj_CdsCt">]
""",
)

entry(
    index = 7,
    label = "CsRHHH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H   u0 {1,S}
3    H   u0 {1,S}
4    H   u0 {1,S}
5    R!H ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.003083, 'd13': 0.009085, 'd23': 0.00437},
        uncertainties = {'d12': 0.159809, 'd13': 0.160865, 'd23': 0.137614},
    ),
    shortDesc = u"""Fitted to 1393 distances.
""",
    longDesc = 
u"""
[<Entry index=12 label="CsOHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=7 label="CsCHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=7 label="CsCHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=172 label="OjO">]
[<Entry index=7 label="CsCHHH">, <Entry index=268 label="CtjC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=7 label="CsCHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=177 label="Cs_methyl">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=268 label="CtjC">]
[<Entry index=7 label="CsCHHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=163 label="Hrad">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=7 label="CsCHHH">, <Entry index=168 label="OjCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=174 label="Crad">]
[<Entry index=7 label="CsCHHH">, <Entry index=172 label="OjO">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=269 label="Cbj">]
[<Entry index=7 label="CsCHHH">, <Entry index=166 label="OjH">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=172 label="OjO">]
[<Entry index=7 label="CsCHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=163 label="Hrad">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=12 label="CsOHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=172 label="OjO">]
[<Entry index=12 label="CsOHHH">, <Entry index=166 label="OjH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=7 label="CsCHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=186 label="CsjCCH">]
[<Entry index=12 label="CsOHHH">, <Entry index=168 label="OjCs">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=174 label="Crad">]
[<Entry index=12 label="CsOHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=166 label="OjH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=269 label="Cbj">]
[<Entry index=7 label="CsCHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=12 label="CsOHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=12 label="CsOHHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=12 label="CsOHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=163 label="Hrad">]
[<Entry index=12 label="CsOHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=179 label="CsjCH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=184 label="CsjOH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=268 label="CtjC">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=204 label="CsjCCC">]
[<Entry index=12 label="CsOHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=12 label="CsOHHH">, <Entry index=269 label="Cbj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=12 label="CsOHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=167 label="OjC">]
[<Entry index=12 label="CsOHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=12 label="CsOHHH">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=12 label="CsOHHH">, <Entry index=163 label="Hrad">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=175 label="Cj">]
[<Entry index=7 label="CsCHHH">, <Entry index=175 label="Cj">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=12 label="CsOHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=7 label="CsCHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=168 label="OjCs">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=172 label="OjO">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=166 label="OjH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=177 label="Cs_methyl">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=7 label="CsCHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=12 label="CsOHHH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=7 label="CsCHHH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=335 label="C_quartetR">]
[<Entry index=7 label="CsCHHH">, <Entry index=269 label="Cbj">]
[<Entry index=12 label="CsOHHH">, <Entry index=181 label="Csj/Cd/H2">]
""",
)

entry(
    index = 8,
    label = "CsCHHH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
5    C  ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.007674, 'd13': 0.016223, 'd23': 0.006095},
        uncertainties = {'d12': 0.160045, 'd13': 0.156463, 'd23': 0.135002},
    ),
    shortDesc = u"""Fitted to 1204 distances.
""",
    longDesc = 
u"""
[<Entry index=8 label="C/H3/Cs">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=7 label="CsCHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=7 label="CsCHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=172 label="OjO">]
[<Entry index=7 label="CsCHHH">, <Entry index=268 label="CtjC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=7 label="CsCHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=268 label="CtjC">]
[<Entry index=7 label="CsCHHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=163 label="Hrad">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=7 label="CsCHHH">, <Entry index=168 label="OjCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=174 label="Crad">]
[<Entry index=7 label="CsCHHH">, <Entry index=172 label="OjO">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=172 label="OjO">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=174 label="Crad">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=269 label="Cbj">]
[<Entry index=7 label="CsCHHH">, <Entry index=166 label="OjH">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=168 label="OjCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=172 label="OjO">]
[<Entry index=7 label="CsCHHH">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=7 label="CsCHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=186 label="CsjCCH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=166 label="OjH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=269 label="Cbj">]
[<Entry index=7 label="CsCHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=163 label="Hrad">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=179 label="CsjCH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=184 label="CsjOH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=268 label="CtjC">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=204 label="CsjCCC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=7 label="CsCHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=175 label="Cj">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=177 label="Cs_methyl">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=7 label="CsCHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=163 label="Hrad">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=166 label="OjH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=177 label="Cs_methyl">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=7 label="CsCHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=7 label="CsCHHH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=335 label="C_quartetR">]
[<Entry index=7 label="CsCHHH">, <Entry index=269 label="Cbj">]
""",
)

entry(
    index = 9,
    label = "C/H3/Cs",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
5    Cs ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.000124, 'd13': 0.000975, 'd23': -0.001576},
        uncertainties = {'d12': 0.16656, 'd13': 0.168268, 'd23': 0.127444},
    ),
    shortDesc = u"""Fitted to 879 distances.
""",
    longDesc = 
u"""
[<Entry index=8 label="C/H3/Cs">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=269 label="Cbj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=268 label="CtjC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=168 label="OjCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=186 label="CsjCCH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=163 label="Hrad">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=179 label="CsjCH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=184 label="CsjOH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=177 label="Cs_methyl">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=204 label="CsjCCC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=166 label="OjH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=172 label="OjO">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=335 label="C_quartetR">]
""",
)

entry(
    index = 10,
    label = "C/H3/Cd",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
5    Cd ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.034554, 'd13': 0.058562, 'd23': 0.02033},
        uncertainties = {'d12': 0.152298, 'd13': 0.111364, 'd23': 0.143687},
    ),
    shortDesc = u"""Fitted to 198 distances.
""",
    longDesc = 
u"""
[<Entry index=9 label="C/H3/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=172 label="OjO">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=166 label="OjH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=269 label="Cbj">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=268 label="CtjC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=175 label="Cj">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=174 label="Crad">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
""",
)

entry(
    index = 11,
    label = "C/H3/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
5    Ct u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.047165, 'd13': 0.063055, 'd23': 0.015511},
        uncertainties = {'d12': 0.149815, 'd13': 0.06748, 'd23': 0.263911},
    ),
    shortDesc = u"""Fitted to 16 distances.
""",
    longDesc = 
u"""
[<Entry index=10 label="C/H3/Ct">, <Entry index=175 label="Cj">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=163 label="Hrad">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=177 label="Cs_methyl">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=172 label="OjO">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=174 label="Crad">]
""",
)

entry(
    index = 12,
    label = "C/H3/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
5    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 13,
    label = "CsOHHH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.024765, 'd13': -0.03421, 'd23': -0.006091},
        uncertainties = {'d12': 0.159753, 'd13': 0.188121, 'd23': 0.154573},
    ),
    shortDesc = u"""Fitted to 189 distances.
""",
    longDesc = 
u"""
[<Entry index=12 label="CsOHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=12 label="CsOHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=167 label="OjC">]
[<Entry index=12 label="CsOHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=12 label="CsOHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=12 label="CsOHHH">, <Entry index=175 label="Cj">]
[<Entry index=12 label="CsOHHH">, <Entry index=163 label="Hrad">]
[<Entry index=12 label="CsOHHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=12 label="CsOHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=12 label="CsOHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=12 label="CsOHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=166 label="OjH">]
[<Entry index=12 label="CsOHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=168 label="OjCs">]
[<Entry index=12 label="CsOHHH">, <Entry index=172 label="OjO">]
[<Entry index=12 label="CsOHHH">, <Entry index=269 label="Cbj">]
[<Entry index=12 label="CsOHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=12 label="CsOHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=12 label="CsOHHH">, <Entry index=250 label="Cdj_CdsCt">]
""",
)

entry(
    index = 14,
    label = "CsRRHH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H   u0 {1,S}
3    H   u0 {1,S}
4    R!H ux {1,S}
5    R!H ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.001934, 'd13': 0.065851, 'd23': 0.061934},
        uncertainties = {'d12': 0.175469, 'd13': 0.102448, 'd23': 0.180196},
    ),
    shortDesc = u"""Fitted to 157 distances.
""",
    longDesc = 
u"""
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=166 label="OjH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=172 label="OjO">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=14 label="CsCCHH">, <Entry index=166 label="OjH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=14 label="CsCCHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=166 label="OjH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=175 label="Cj">]
[<Entry index=14 label="CsCCHH">, <Entry index=163 label="Hrad">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=172 label="OjO">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=166 label="OjH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=163 label="Hrad">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=175 label="Cj">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=172 label="OjO">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=204 label="CsjCCC">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=166 label="OjH">]
[<Entry index=14 label="CsCCHH">, <Entry index=172 label="OjO">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=14 label="CsCCHH">, <Entry index=168 label="OjCs">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=172 label="OjO">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=168 label="OjCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=25 label="CsCOHH">, <Entry index=163 label="Hrad">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=163 label="Hrad">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=25 label="CsCOHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=15 label="C/H2/Cs/Cs">, <Entry index=172 label="OjO">]
[<Entry index=25 label="CsCOHH">, <Entry index=172 label="OjO">]
""",
)

entry(
    index = 15,
    label = "CsCCHH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    C  ux {1,S}
5    C  ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.001475, 'd13': 0.074974, 'd23': 0.07144},
        uncertainties = {'d12': 0.171819, 'd13': 0.075615, 'd23': 0.18538},
    ),
    shortDesc = u"""Fitted to 136 distances.
""",
    longDesc = 
u"""
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=172 label="OjO">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=14 label="CsCCHH">, <Entry index=166 label="OjH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=166 label="OjH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=175 label="Cj">]
[<Entry index=14 label="CsCCHH">, <Entry index=163 label="Hrad">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=172 label="OjO">]
[<Entry index=15 label="C/H2/Cs/Cs">, <Entry index=172 label="OjO">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=175 label="Cj">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=14 label="CsCCHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=204 label="CsjCCC">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=166 label="OjH">]
[<Entry index=14 label="CsCCHH">, <Entry index=172 label="OjO">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=14 label="CsCCHH">, <Entry index=168 label="OjCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=300 label="Cs_tripH2">]
""",
)

entry(
    index = 16,
    label = "C/H2/Cs/Cs",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cs ux {1,S}
5    Cs ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.091974, 'd13': -0.113981, 'd23': -0.038984},
        uncertainties = {},
    ),
    shortDesc = u"""Fitted to 2 distances.
""",
    longDesc = 
u"""
[<Entry index=15 label="C/H2/Cs/Cs">, <Entry index=172 label="OjO">]
""",
)

entry(
    index = 17,
    label = "C/H2/Cs/Cd",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cs ux {1,S}
5    Cd ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.022394, 'd13': 0.078091, 'd23': 0.100302},
        uncertainties = {'d12': 0.089471, 'd13': 0.048721, 'd23': 0.096706},
    ),
    shortDesc = u"""Fitted to 44 distances.
""",
    longDesc = 
u"""
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=204 label="CsjCCC">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=175 label="Cj">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=166 label="OjH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=172 label="OjO">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=181 label="Csj/Cd/H2">]
""",
)

entry(
    index = 18,
    label = "C/H2/Cs/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cs ux {1,S}
5    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 19,
    label = "C/H2/Cs/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cs ux {1,S}
5    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 20,
    label = "C/H2/Cd/Cd",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cd ux {1,S}
5    Cd ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.015707, 'd13': 0.076655, 'd23': 0.05804},
        uncertainties = {'d12': 0.207913, 'd13': 0.089383, 'd23': 0.22473},
    ),
    shortDesc = u"""Fitted to 83 distances.
""",
    longDesc = 
u"""
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=175 label="Cj">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=172 label="OjO">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=166 label="OjH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=163 label="Hrad">]
""",
)

entry(
    index = 21,
    label = "C/H2/Cd/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cd ux {1,S}
5    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 22,
    label = "C/H2/Cd/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cd ux {1,S}
5    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 23,
    label = "C/H2/Ct/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Ct u0 {1,S}
5    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 24,
    label = "C/H2/Ct/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Ct u0 {1,S}
5    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 25,
    label = "C/H2/Cb/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cb u0 {1,S}
5    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 26,
    label = "CsCOHH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    C  ux {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.005609, 'd13': -0.007133, 'd23': -0.01411},
        uncertainties = {'d12': 0.214788, 'd13': 0.219914, 'd23': 0.156158},
    ),
    shortDesc = u"""Fitted to 21 distances.
""",
    longDesc = 
u"""
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=166 label="OjH">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=172 label="OjO">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=166 label="OjH">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=163 label="Hrad">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=168 label="OjCs">]
[<Entry index=25 label="CsCOHH">, <Entry index=163 label="Hrad">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=172 label="OjO">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=163 label="Hrad">]
[<Entry index=25 label="CsCOHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=25 label="CsCOHH">, <Entry index=172 label="OjO">]
""",
)

entry(
    index = 27,
    label = "C/H2/Cs/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cs ux {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.029363, 'd13': -0.012396, 'd23': -0.040152},
        uncertainties = {'d12': 0.078209, 'd13': 0.087225, 'd23': 0.022803},
    ),
    shortDesc = u"""Fitted to 7 distances.
""",
    longDesc = 
u"""
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=163 label="Hrad">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=166 label="OjH">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=172 label="OjO">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=168 label="OjCs">]
""",
)

entry(
    index = 28,
    label = "C/H2/Cd/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cd ux {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.013225, 'd13': -0.033503, 'd23': -0.022173},
        uncertainties = {'d12': 0.33764, 'd13': 0.325142, 'd23': 0.152692},
    ),
    shortDesc = u"""Fitted to 9 distances.
""",
    longDesc = 
u"""
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=172 label="OjO">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=166 label="OjH">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=163 label="Hrad">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=177 label="Cs_methyl">]
""",
)

entry(
    index = 29,
    label = "C/H2/Ct/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Ct u0 {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 30,
    label = "C/H2/Cb/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cb u0 {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 31,
    label = "CsOOHH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    O  ux {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
    shortDesc = u"""Fitted to 2 distances.
""",
)

entry(
    index = 32,
    label = "CsRRRH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H   u0 {1,S}
3    R!H ux {1,S}
4    R!H ux {1,S}
5    R!H ux {1,S}
""",
    kinetics = DistanceData(distances={}),
    shortDesc = u"""Fitted to 37 distances.
""",
)

entry(
    index = 33,
    label = "CsCCCH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    C  ux {1,S}
4    C  ux {1,S}
5    C  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
    shortDesc = u"""Fitted to 30 distances.
""",
)

entry(
    index = 34,
    label = "C/H/Cs/Cs/Cs",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    Cs ux {1,S}
5    Cs ux {1,S}
""",
    kinetics = DistanceData(distances={}),
    shortDesc = u"""Fitted to 13 distances.
""",
)

entry(
    index = 35,
    label = "C/H/Cs/Cs/Cd",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    Cs ux {1,S}
5    Cd ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 36,
    label = "C/H/Cs/Cs/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    Cs ux {1,S}
5    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 37,
    label = "C/H/Cs/Cs/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    Cs ux {1,S}
5    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 38,
    label = "C/H/Cs/Cd/Cd",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    Cd ux {1,S}
5    Cd ux {1,S}
""",
    kinetics = DistanceData(distances={}),
    shortDesc = u"""Fitted to 1 distances.
""",
)

entry(
    index = 39,
    label = "C/H/Cs/Cd/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    Cd ux {1,S}
5    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 40,
    label = "C/H/Cs/Cd/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    Cd ux {1,S}
5    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 41,
    label = "C/H/Cs/Ct/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    Ct u0 {1,S}
5    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 42,
    label = "C/H/Cs/Ct/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    Ct u0 {1,S}
5    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 43,
    label = "C/H/Cs/Cb/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    Cb u0 {1,S}
5    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 44,
    label = "C/H/Cd/Cd/Cd",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,S}
4    Cd ux {1,S}
5    Cd ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 45,
    label = "C/H/Cd/Cd/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,S}
4    Cd ux {1,S}
5    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 46,
    label = "C/H/Cd/Cd/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,S}
4    Cd ux {1,S}
5    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 47,
    label = "C/H/Cd/Ct/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,S}
4    Ct u0 {1,S}
5    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 48,
    label = "C/H/Cd/Ct/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,S}
4    Ct u0 {1,S}
5    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 49,
    label = "C/H/Cd/Cb/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,S}
4    Cb u0 {1,S}
5    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 50,
    label = "C/H/Ct/Ct/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Ct u0 {1,S}
4    Ct u0 {1,S}
5    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 51,
    label = "C/H/Ct/Ct/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Ct u0 {1,S}
4    Ct u0 {1,S}
5    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 52,
    label = "C/H/Ct/Cb/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Ct u0 {1,S}
4    Cb u0 {1,S}
5    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 53,
    label = "C/H/Cb/Cb/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cb u0 {1,S}
4    Cb u0 {1,S}
5    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 54,
    label = "CsCCOH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    C  ux {1,S}
4    C  ux {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
    shortDesc = u"""Fitted to 7 distances.
""",
)

entry(
    index = 55,
    label = "C/H/Cs/Cs/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    Cs ux {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
    shortDesc = u"""Fitted to 7 distances.
""",
)

entry(
    index = 56,
    label = "C/H/Cs/Cd/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    Cd ux {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 57,
    label = "C/H/Cs/Ct/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    Ct u0 {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 58,
    label = "C/H/Cs/Cb/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    Cb u0 {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 59,
    label = "C/H/Cd/Cd/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,S}
4    Cd ux {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 60,
    label = "C/H/Cd/Ct/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,S}
4    Ct u0 {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 61,
    label = "C/H/Cd/Cb/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,S}
4    Cb u0 {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 62,
    label = "C/H/Ct/Ct/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Ct u0 {1,S}
4    Ct u0 {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 63,
    label = "C/H/Ct/Cb/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Ct u0 {1,S}
4    Cb u0 {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 64,
    label = "C/H/Cb/Cb/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cb u0 {1,S}
4    Cb u0 {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 65,
    label = "CsCOOH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    C  ux {1,S}
4    O  ux {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 66,
    label = "C/H/Cs/O/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    O  ux {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 67,
    label = "C/H/Cd/O/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,S}
4    O  ux {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 68,
    label = "C/H/Ct/O/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Ct u0 {1,S}
4    O  ux {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 69,
    label = "C/H/Cb/O/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cb u0 {1,S}
4    O  ux {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 70,
    label = "CsOOOH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    O  ux {1,S}
4    O  ux {1,S}
5    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 71,
    label = "Csrad_H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S}
2 *2 H  u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.068378, 'd13': 0.040085, 'd23': -0.03457},
        uncertainties = {'d12': 0.108092, 'd13': 0.073317, 'd23': 0.058283},
    ),
    shortDesc = u"""Fitted to 27 distances.
""",
    longDesc = 
u"""
[<Entry index=71 label="C_methyl">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=175 label="Cj">]
[<Entry index=71 label="C_methyl">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=163 label="Hrad">]
[<Entry index=71 label="C_methyl">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=71 label="C_methyl">, <Entry index=204 label="CsjCCC">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=269 label="Cbj">]
[<Entry index=71 label="C_methyl">, <Entry index=172 label="OjO">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=163 label="Hrad">]
[<Entry index=71 label="C_methyl">, <Entry index=168 label="OjCs">]
[<Entry index=71 label="C_methyl">, <Entry index=166 label="OjH">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=71 label="C_methyl">, <Entry index=175 label="Cj">]
[<Entry index=71 label="C_methyl">, <Entry index=186 label="CsjCCH">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=166 label="OjH">]
[<Entry index=71 label="C_methyl">, <Entry index=188 label="Csj/Cs/Cd/H">]
""",
)

entry(
    index = 72,
    label = "C_methyl",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.097011, 'd13': 0.04669, 'd23': -0.057816},
        uncertainties = {'d12': 0.126678, 'd13': 0.091589, 'd23': 0.05457},
    ),
    shortDesc = u"""Fitted to 17 distances.
""",
    longDesc = 
u"""
[<Entry index=71 label="C_methyl">, <Entry index=168 label="OjCs">]
[<Entry index=71 label="C_methyl">, <Entry index=166 label="OjH">]
[<Entry index=71 label="C_methyl">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=71 label="C_methyl">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=71 label="C_methyl">, <Entry index=186 label="CsjCCH">]
[<Entry index=71 label="C_methyl">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=163 label="Hrad">]
[<Entry index=71 label="C_methyl">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=204 label="CsjCCC">]
[<Entry index=71 label="C_methyl">, <Entry index=172 label="OjO">]
[<Entry index=71 label="C_methyl">, <Entry index=175 label="Cj">]
""",
)

entry(
    index = 73,
    label = "CsradRH2",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs  u1 {2,S} {3,S} {4,S}
2 *2 H   u0 {1,S}
3    H   u0 {1,S}
4    R!H ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.008726, 'd13': 0.026324, 'd23': 0.013858},
        uncertainties = {'d12': 0.091622, 'd13': 0.042736, 'd23': 0.076531},
    ),
    shortDesc = u"""Fitted to 10 distances.
""",
    longDesc = 
u"""
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=163 label="Hrad">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=166 label="OjH">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=175 label="Cj">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=269 label="Cbj">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=179 label="CsjCH2">]
""",
)

entry(
    index = 74,
    label = "CsradCHH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    C  ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.008726, 'd13': 0.026324, 'd23': 0.013858},
        uncertainties = {'d12': 0.091622, 'd13': 0.042736, 'd23': 0.076531},
    ),
    shortDesc = u"""Fitted to 10 distances.
""",
    longDesc = 
u"""
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=163 label="Hrad">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=166 label="OjH">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=175 label="Cj">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=269 label="Cbj">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=179 label="CsjCH2">]
""",
)

entry(
    index = 75,
    label = "Csrad/H/Cs/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cs ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 76,
    label = "Csrad/H/Cd/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cd ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 77,
    label = "Csrad/H/Ct/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Ct u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.008726, 'd13': 0.026324, 'd23': 0.013858},
        uncertainties = {'d12': 0.091622, 'd13': 0.042736, 'd23': 0.076531},
    ),
    shortDesc = u"""Fitted to 10 distances.
""",
    longDesc = 
u"""
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=163 label="Hrad">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=166 label="OjH">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=175 label="Cj">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=269 label="Cbj">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=179 label="CsjCH2">]
""",
)

entry(
    index = 78,
    label = "Csrad/H/Cb/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 79,
    label = "CsradOH2",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 80,
    label = "CsradRRH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs  u1 {2,S} {3,S} {4,S}
2 *2 H   u0 {1,S}
3    R!H ux {1,S}
4    R!H ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 81,
    label = "CsradCCH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    C  ux {1,S}
4    C  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 82,
    label = "Csrad/Cs/Cs/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    Cs ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 83,
    label = "Csrad/Cs/Cd/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    Cd ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 84,
    label = "Csrad/Cs/Ct/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 85,
    label = "Csrad/Cs/Cb/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 86,
    label = "Csrad/Cd/Cd/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,S}
4    Cd ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 87,
    label = "Csrad/Cd/Ct/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,S}
4    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 88,
    label = "Csrad/Cd/Cb/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 89,
    label = "Csrad/Ct/Ct/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    Ct u0 {1,S}
4    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 90,
    label = "Csrad/Ct/Cb/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    Ct u0 {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 91,
    label = "Csrad/Cb/Cb/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    Cb u0 {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 92,
    label = "CsradCOH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    C  ux {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 93,
    label = "Csrad/Cs/O/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 94,
    label = "Csrad/Cd/O/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 95,
    label = "Csrad/Ct/O/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    Ct u0 {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 96,
    label = "Csrad/Cb/O/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    Cb u0 {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 97,
    label = "CsradOOH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    O  ux {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 98,
    label = "CsbiradH",
    group = "OR{Cs_singletH, Cs_tripletH}",
    kinetics = DistanceData(
        distances = {'d12': 0.081422, 'd13': -0.069261, 'd23': -0.13074},
        uncertainties = {'d12': 0.357824, 'd13': 0.435552, 'd23': 0.236232},
    ),
    shortDesc = u"""Fitted to 3 distances.
""",
    longDesc = 
u"""
[<Entry index=108 label="Cs_tripletHH">, <Entry index=163 label="Hrad">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=166 label="OjH">]
""",
)

entry(
    index = 99,
    label = "Cs_singletH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 p1 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    R  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 100,
    label = "Cs_singletHH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 p1 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 101,
    label = "Cs_singletRH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs  u0 p1 {2,S} {3,S}
2 *2 H   u0 {1,S}
3    R!H ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 102,
    label = "C_singletCH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 p1 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    C  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 103,
    label = "C_singlet/Cs/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 p1 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 104,
    label = "C_singlet/Cd/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 p1 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 105,
    label = "C_singlet/Ct/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 p1 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 106,
    label = "C_singlet/Cb/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 p1 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 107,
    label = "C_singletOH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u0 p1 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 108,
    label = "Cs_tripletH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u2 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    R  ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.081422, 'd13': -0.069261, 'd23': -0.13074},
        uncertainties = {'d12': 0.357824, 'd13': 0.435552, 'd23': 0.236232},
    ),
    shortDesc = u"""Fitted to 3 distances.
""",
    longDesc = 
u"""
[<Entry index=108 label="Cs_tripletHH">, <Entry index=163 label="Hrad">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=166 label="OjH">]
""",
)

entry(
    index = 109,
    label = "Cs_tripletHH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u2 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.081422, 'd13': -0.069261, 'd23': -0.13074},
        uncertainties = {'d12': 0.357824, 'd13': 0.435552, 'd23': 0.236232},
    ),
    shortDesc = u"""Fitted to 3 distances.
""",
    longDesc = 
u"""
[<Entry index=108 label="Cs_tripletHH">, <Entry index=163 label="Hrad">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=166 label="OjH">]
""",
)

entry(
    index = 110,
    label = "Cs_tripletRH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs  u2 {2,S} {3,S}
2 *2 H   u0 {1,S}
3    R!H ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 111,
    label = "Cs_tripletCH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u2 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    C  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 112,
    label = "C_triplet/Cs/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u2 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 113,
    label = "C_triplet/Cd/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u2 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 114,
    label = "C_triplet/Ct/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u2 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 115,
    label = "C_triplet/Cb/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u2 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 116,
    label = "Cs_tripletOH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cs u2 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 117,
    label = "CstriradH",
    group = "OR{Cdoublet_H, Cquartet_H}",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 118,
    label = "Cdoublet_H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 C u1 p1 {2,S}
2 *2 H u0 p0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 119,
    label = "Cquartet_H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 C u3 p0 {2,S}
2 *2 H u0 p0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 120,
    label = "Cd_H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd ux {2,S}
2 *2 H  u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.064316, 'd13': 0.048299, 'd23': -0.016606},
        uncertainties = {'d12': 0.152976, 'd13': 0.11197, 'd23': 0.13496},
    ),
    shortDesc = u"""Fitted to 270 distances.
""",
    longDesc = 
u"""
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=172 label="OjO">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=269 label="Cbj">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=168 label="OjCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=166 label="OjH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=175 label="Cj">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=186 label="CsjCCH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=172 label="OjO">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=186 label="CsjCCH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=204 label="CsjCCC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=174 label="Crad">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=127 label="Cd_Cds/Cd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=163 label="Hrad">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=204 label="CsjCCC">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=269 label="Cbj">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=168 label="OjCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=166 label="OjH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=175 label="Cj">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=184 label="CsjOH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=172 label="OjO">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=268 label="CtjC">]
[<Entry index=121 label="Cd_C/R/H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=163 label="Hrad">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
""",
)

entry(
    index = 121,
    label = "Cdnorad_H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u0 {2,S}
2 *2 H  u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.068626, 'd13': 0.050909, 'd23': -0.019505},
        uncertainties = {'d12': 0.15776, 'd13': 0.115234, 'd23': 0.139255},
    ),
    shortDesc = u"""Fitted to 248 distances.
""",
    longDesc = 
u"""
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=269 label="Cbj">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=168 label="OjCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=166 label="OjH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=175 label="Cj">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=186 label="CsjCCH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=172 label="OjO">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=204 label="CsjCCC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=174 label="Crad">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=127 label="Cd_Cds/Cd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=163 label="Hrad">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=269 label="Cbj">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=168 label="OjCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=166 label="OjH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=175 label="Cj">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=172 label="OjO">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=268 label="CtjC">]
[<Entry index=121 label="Cd_C/R/H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=163 label="Hrad">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
""",
)

entry(
    index = 122,
    label = "Cd_C/R/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u0 {2,S} {3,D}
2 *2 H  u0 {1,S}
3    C  ux {1,D}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.07427, 'd13': 0.050886, 'd23': -0.024913},
        uncertainties = {'d12': 0.159538, 'd13': 0.114794, 'd23': 0.139928},
    ),
    shortDesc = u"""Fitted to 236 distances.
""",
    longDesc = 
u"""
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=269 label="Cbj">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=168 label="OjCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=166 label="OjH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=175 label="Cj">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=186 label="CsjCCH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=172 label="OjO">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=204 label="CsjCCC">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=174 label="Crad">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=127 label="Cd_Cds/Cd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=163 label="Hrad">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=269 label="Cbj">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=168 label="OjCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=166 label="OjH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=175 label="Cj">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=163 label="Hrad">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=172 label="OjO">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=268 label="CtjC">]
[<Entry index=121 label="Cd_C/R/H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
""",
)

entry(
    index = 123,
    label = "Cd_C/H2",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u0 {2,S} {3,D} {4,S}
2 *2 H  u0 {1,S}
3    C  ux {1,D}
4    H  u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.077902, 'd13': 0.050732, 'd23': -0.028754},
        uncertainties = {'d12': 0.153779, 'd13': 0.114541, 'd23': 0.136592},
    ),
    shortDesc = u"""Fitted to 229 distances.
""",
    longDesc = 
u"""
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=269 label="Cbj">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=168 label="OjCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=166 label="OjH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=175 label="Cj">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=186 label="CsjCCH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=172 label="OjO">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=204 label="CsjCCC">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=174 label="Crad">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=163 label="Hrad">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=269 label="Cbj">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=168 label="OjCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=166 label="OjH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=175 label="Cj">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=163 label="Hrad">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=172 label="OjO">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=268 label="CtjC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
""",
)

entry(
    index = 124,
    label = "Cd_Cds/H2",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u0 {2,S} {3,D} {4,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,D}
4    H  u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.076501, 'd13': 0.054188, 'd23': -0.024801},
        uncertainties = {'d12': 0.163976, 'd13': 0.122219, 'd23': 0.147778},
    ),
    shortDesc = u"""Fitted to 172 distances.
""",
    longDesc = 
u"""
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=168 label="OjCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=166 label="OjH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=186 label="CsjCCH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=174 label="Crad">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=163 label="Hrad">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=269 label="Cbj">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=175 label="Cj">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=172 label="OjO">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=268 label="CtjC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=246 label="Cdj_CddH">]
""",
)

entry(
    index = 125,
    label = "Cd_Cdd/H2",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd  u0 {2,S} {3,D} {4,S}
2 *2 H   u0 {1,S}
3    Cdd u0 {1,D}
4    H   u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.081965, 'd13': 0.040713, 'd23': -0.040211},
        uncertainties = {'d12': 0.122283, 'd13': 0.090733, 'd23': 0.099166},
    ),
    shortDesc = u"""Fitted to 57 distances.
""",
    longDesc = 
u"""
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=269 label="Cbj">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=168 label="OjCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=166 label="OjH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=175 label="Cj">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=204 label="CsjCCC">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=172 label="OjO">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=163 label="Hrad">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=331 label="Cd_tripletC">]
""",
)

entry(
    index = 126,
    label = "Cd_C/C/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u0 {2,S} {3,D} {4,S}
2 *2 H  u0 {1,S}
3    C  ux {1,D}
4    C  ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.073946, 'd13': 0.043002, 'd23': 0.110128},
        uncertainties = {},
    ),
    shortDesc = u"""Fitted to 2 distances.
""",
    longDesc = 
u"""
[<Entry index=127 label="Cd_Cds/Cd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
""",
)

entry(
    index = 127,
    label = "Cd_Cds/Cs/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u0 {2,S} {3,D} {4,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,D}
4    Cs ux {1,S}
""",
    kinetics = DistanceData(distances={}),
    shortDesc = u"""Fitted to 40 distances.
""",
)

entry(
    index = 128,
    label = "Cd_Cds/Cd/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u0 {2,S} {3,D} {4,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,D}
4    Cd ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.073946, 'd13': 0.043002, 'd23': 0.110128},
        uncertainties = {},
    ),
    shortDesc = u"""Fitted to 2 distances.
""",
    longDesc = 
u"""
[<Entry index=127 label="Cd_Cds/Cd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
""",
)

entry(
    index = 129,
    label = "Cd_Cds/Ct/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u0 {2,S} {3,D} {4,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,D}
4    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
    shortDesc = u"""Fitted to 15 distances.
""",
)

entry(
    index = 130,
    label = "Cd_Cds/Cb/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u0 {2,S} {3,D} {4,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,D}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 131,
    label = "Cd_Cdd/Cs/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd  u0 {2,S} {3,D} {4,S}
2 *2 H   u0 {1,S}
3    Cdd u0 {1,D}
4    Cs  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 132,
    label = "Cd_Cdd/Cd/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd  u0 {2,S} {3,D} {4,S}
2 *2 H   u0 {1,S}
3    Cdd u0 {1,D}
4    Cd  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 133,
    label = "Cd_Cdd/Ct/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd  u0 {2,S} {3,D} {4,S}
2 *2 H   u0 {1,S}
3    Cdd u0 {1,D}
4    Ct  u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 134,
    label = "Cd_Cdd/Cb/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd  u0 {2,S} {3,D} {4,S}
2 *2 H   u0 {1,S}
3    Cdd u0 {1,D}
4    Cb  u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 135,
    label = "Cd_C/O/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u0 {2,S} {3,D} {4,S}
2 *2 H  u0 {1,S}
3    C  ux {1,D}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 136,
    label = "Cd_Cds/O/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u0 {2,S} {3,D} {4,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,D}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 137,
    label = "Cd_Cdd/O/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd  u0 {2,S} {3,D} {4,S}
2 *2 H   u0 {1,S}
3    Cdd ux {1,D}
4    O   ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 138,
    label = "Cd_O/R/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u0 {2,S} {3,D}
2 *2 H  u0 {1,S}
3    O  u0 {1,D}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 139,
    label = "Cd_O/H2",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u0 {2,S} {3,D} {4,S}
2 *2 H  u0 {1,S}
3    O  u0 {1,D}
4    H  u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 140,
    label = "Cd_O/C/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u0 {2,S} {3,D} {4,S}
2 *2 H  u0 {1,S}
3    O  u0 {1,D}
4    C  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 141,
    label = "Cd_O/Cs/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u0 {2,S} {3,D} {4,S}
2 *2 H  u0 {1,S}
3    O  u0 {1,D}
4    Cs ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 142,
    label = "Cd_O/Cd/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u0 {2,S} {3,D} {4,S}
2 *2 H  u0 {1,S}
3    O  u0 {1,D}
4    Cd ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 143,
    label = "Cd_O/Ct/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u0 {2,S} {3,D} {4,S}
2 *2 H  u0 {1,S}
3    O  u0 {1,D}
4    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 144,
    label = "Cd_O/Cb/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u0 {2,S} {3,D} {4,S}
2 *2 H  u0 {1,S}
3    O  u0 {1,D}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 145,
    label = "Cd_O/O/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u0 {2,S} {3,D} {4,S}
2 *2 H  u0 {1,S}
3    O  u0 {1,D}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 146,
    label = "Cdrad_H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u1 {2,S}
2 *2 H  u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.018467, 'd13': 0.020529, 'd23': 0.014242},
        uncertainties = {'d12': 0.090045, 'd13': 0.071138, 'd23': 0.077696},
    ),
    shortDesc = u"""Fitted to 22 distances.
""",
    longDesc = 
u"""
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=184 label="CsjOH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=172 label="OjO">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=186 label="CsjCCH">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=204 label="CsjCCC">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=180 label="Csj/Cs/H2">]
""",
)

entry(
    index = 147,
    label = "Cdrad_C/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u1 {2,S} {3,D}
2 *2 H  u0 {1,S}
3    C  ux {1,D}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.018467, 'd13': 0.020529, 'd23': 0.014242},
        uncertainties = {'d12': 0.090045, 'd13': 0.071138, 'd23': 0.077696},
    ),
    shortDesc = u"""Fitted to 22 distances.
""",
    longDesc = 
u"""
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=184 label="CsjOH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=172 label="OjO">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=186 label="CsjCCH">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=204 label="CsjCCC">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=180 label="Csj/Cs/H2">]
""",
)

entry(
    index = 148,
    label = "Cdrad_Cds/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u1 {2,S} {3,D}
2 *2 H  u0 {1,S}
3    Cd ux {1,D}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 149,
    label = "Cdrad_Cdd/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd  u1 {2,S} {3,D}
2 *2 H   u0 {1,S}
3    Cdd u0 {1,D}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.018467, 'd13': 0.020529, 'd23': 0.014242},
        uncertainties = {'d12': 0.090045, 'd13': 0.071138, 'd23': 0.077696},
    ),
    shortDesc = u"""Fitted to 22 distances.
""",
    longDesc = 
u"""
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=184 label="CsjOH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=172 label="OjO">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=186 label="CsjCCH">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=204 label="CsjCCC">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=180 label="Csj/Cs/H2">]
""",
)

entry(
    index = 150,
    label = "Cdrad_O/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cd u1 {2,S} {3,D}
2 *2 H  u0 {1,S}
3    O  u0 {1,D}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 151,
    label = "Ct_H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Ct u0 {2,S}
2 *2 H  u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.25186, 'd13': 0.174208, 'd23': -0.027581},
        uncertainties = {'d12': 0.437512, 'd13': 0.245291, 'd23': 0.311807},
    ),
    shortDesc = u"""Fitted to 39 distances.
""",
    longDesc = 
u"""
[<Entry index=150 label="Ct_H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=150 label="Ct_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=150 label="Ct_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=150 label="Ct_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=150 label="Ct_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=150 label="Ct_H">, <Entry index=268 label="CtjC">]
[<Entry index=150 label="Ct_H">, <Entry index=269 label="Cbj">]
[<Entry index=150 label="Ct_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=150 label="Ct_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=150 label="Ct_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
""",
)

entry(
    index = 152,
    label = "Cb_H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 Cb u0 {2,S}
2 *2 H  u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.077476, 'd13': 0.042539, 'd23': -0.039775},
        uncertainties = {'d12': 0.234691, 'd13': 0.156054, 'd23': 0.164275},
    ),
    shortDesc = u"""Fitted to 44 distances.
""",
    longDesc = 
u"""
[<Entry index=151 label="Cb_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=151 label="Cb_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=151 label="Cb_H">, <Entry index=166 label="OjH">]
[<Entry index=151 label="Cb_H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=151 label="Cb_H">, <Entry index=175 label="Cj">]
[<Entry index=151 label="Cb_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=151 label="Cb_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=172 label="OjO">]
[<Entry index=151 label="Cb_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=151 label="Cb_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=151 label="Cb_H">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=151 label="Cb_H">, <Entry index=269 label="Cbj">]
[<Entry index=151 label="Cb_H">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=151 label="Cb_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=163 label="Hrad">]
[<Entry index=151 label="Cb_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=151 label="Cb_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=151 label="Cb_H">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=151 label="Cb_H">, <Entry index=184 label="CsjOH2">]
""",
)

entry(
    index = 153,
    label = "O_H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 O ux {2,S}
2 *2 H u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.035947, 'd13': -0.072809, 'd23': -0.02967},
        uncertainties = {'d12': 0.141589, 'd13': 0.134551, 'd23': 0.165049},
    ),
    shortDesc = u"""Fitted to 243 distances.
""",
    longDesc = 
u"""
[<Entry index=153 label="OradH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=161 label="OOH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=161 label="OOH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=155 label="OHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=161 label="OOH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=153 label="OradH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=153 label="OradH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=155 label="OHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=155 label="OHH">, <Entry index=269 label="Cbj">]
[<Entry index=155 label="OHH">, <Entry index=168 label="OjCs">]
[<Entry index=161 label="OOH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=153 label="OradH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=175 label="Cj">]
[<Entry index=153 label="OradH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=153 label="OradH">, <Entry index=179 label="CsjCH2">]
[<Entry index=161 label="OOH">, <Entry index=186 label="CsjCCH">]
[<Entry index=153 label="OradH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=155 label="OHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=155 label="OHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=153 label="OradH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=161 label="OOH">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=155 label="OHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=155 label="OHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=155 label="OHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=153 label="OradH">, <Entry index=204 label="CsjCCC">]
[<Entry index=161 label="OOH">, <Entry index=184 label="CsjOH2">]
[<Entry index=153 label="OradH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=155 label="OHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=155 label="OHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=161 label="OOH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=155 label="OHH">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=153 label="OradH">, <Entry index=184 label="CsjOH2">]
[<Entry index=153 label="OradH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=155 label="OHH">, <Entry index=175 label="Cj">]
[<Entry index=161 label="OOH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=155 label="OHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=161 label="OOH">, <Entry index=169 label="OjCd">]
[<Entry index=154 label="ORH">, <Entry index=184 label="CsjOH2">]
[<Entry index=155 label="OHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=163 label="Hrad">]
[<Entry index=153 label="OradH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=161 label="OOH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=153 label="OradH">, <Entry index=172 label="OjO">]
[<Entry index=153 label="OradH">, <Entry index=168 label="OjCs">]
[<Entry index=161 label="OOH">, <Entry index=175 label="Cj">]
[<Entry index=161 label="OOH">, <Entry index=172 label="OjO">]
[<Entry index=155 label="OHH">, <Entry index=335 label="C_quartetR">]
[<Entry index=161 label="OOH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=155 label="OHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=154 label="ORH">, <Entry index=172 label="OjO">]
[<Entry index=161 label="OOH">, <Entry index=163 label="Hrad">]
[<Entry index=161 label="OOH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=155 label="OHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=161 label="OOH">, <Entry index=179 label="CsjCH2">]
[<Entry index=155 label="OHH">, <Entry index=172 label="OjO">]
[<Entry index=153 label="OradH">, <Entry index=163 label="Hrad">]
[<Entry index=153 label="OradH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=161 label="OOH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=153 label="OradH">, <Entry index=186 label="CsjCCH">]
[<Entry index=161 label="OOH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=155 label="OHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=168 label="OjCs">]
[<Entry index=155 label="OHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=161 label="OOH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=153 label="OradH">, <Entry index=191 label="Csj/Cd/Cd/H">]
""",
)

entry(
    index = 154,
    label = "OradH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 O u1 {2,S}
2 *2 H u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.032793, 'd13': -0.110713, 'd23': -0.082314},
        uncertainties = {'d12': 0.147906, 'd13': 0.099534, 'd23': 0.076879},
    ),
    shortDesc = u"""Fitted to 42 distances.
""",
    longDesc = 
u"""
[<Entry index=153 label="OradH">, <Entry index=179 label="CsjCH2">]
[<Entry index=153 label="OradH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=153 label="OradH">, <Entry index=186 label="CsjCCH">]
[<Entry index=153 label="OradH">, <Entry index=184 label="CsjOH2">]
[<Entry index=153 label="OradH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=153 label="OradH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=153 label="OradH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=153 label="OradH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=153 label="OradH">, <Entry index=204 label="CsjCCC">]
[<Entry index=153 label="OradH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=153 label="OradH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=153 label="OradH">, <Entry index=163 label="Hrad">]
[<Entry index=153 label="OradH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=153 label="OradH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=153 label="OradH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=153 label="OradH">, <Entry index=172 label="OjO">]
[<Entry index=153 label="OradH">, <Entry index=168 label="OjCs">]
[<Entry index=153 label="OradH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=153 label="OradH">, <Entry index=175 label="Cj">]
[<Entry index=153 label="OradH">, <Entry index=191 label="Csj/Cd/Cd/H">]
""",
)

entry(
    index = 155,
    label = "ORH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 O u0 {2,S}
2 *2 H u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.036607, 'd13': -0.064883, 'd23': -0.018661},
        uncertainties = {'d12': 0.141538, 'd13': 0.141478, 'd23': 0.178539},
    ),
    shortDesc = u"""Fitted to 201 distances.
""",
    longDesc = 
u"""
[<Entry index=161 label="OOH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=161 label="OOH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=155 label="OHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=161 label="OOH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=155 label="OHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=155 label="OHH">, <Entry index=269 label="Cbj">]
[<Entry index=155 label="OHH">, <Entry index=168 label="OjCs">]
[<Entry index=161 label="OOH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=161 label="OOH">, <Entry index=184 label="CsjOH2">]
[<Entry index=155 label="OHH">, <Entry index=163 label="Hrad">]
[<Entry index=161 label="OOH">, <Entry index=186 label="CsjCCH">]
[<Entry index=155 label="OHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=155 label="OHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=161 label="OOH">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=155 label="OHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=155 label="OHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=155 label="OHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=155 label="OHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=155 label="OHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=161 label="OOH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=155 label="OHH">, <Entry index=175 label="Cj">]
[<Entry index=155 label="OHH">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=155 label="OHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=161 label="OOH">, <Entry index=169 label="OjCd">]
[<Entry index=154 label="ORH">, <Entry index=184 label="CsjOH2">]
[<Entry index=155 label="OHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=161 label="OOH">, <Entry index=175 label="Cj">]
[<Entry index=161 label="OOH">, <Entry index=172 label="OjO">]
[<Entry index=155 label="OHH">, <Entry index=335 label="C_quartetR">]
[<Entry index=161 label="OOH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=155 label="OHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=154 label="ORH">, <Entry index=172 label="OjO">]
[<Entry index=161 label="OOH">, <Entry index=163 label="Hrad">]
[<Entry index=161 label="OOH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=155 label="OHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=161 label="OOH">, <Entry index=179 label="CsjCH2">]
[<Entry index=155 label="OHH">, <Entry index=172 label="OjO">]
[<Entry index=161 label="OOH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=161 label="OOH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=155 label="OHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=161 label="OOH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=161 label="OOH">, <Entry index=168 label="OjCs">]
""",
)

entry(
    index = 156,
    label = "OHH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 O u0 {2,S} {3,S}
2 *2 H u0 {1,S}
3    H u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.124178, 'd13': -0.063674, 'd23': -0.167425},
        uncertainties = {'d12': 0.159269, 'd13': 0.159608, 'd23': 0.068749},
    ),
    shortDesc = u"""Fitted to 74 distances.
""",
    longDesc = 
u"""
[<Entry index=155 label="OHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=155 label="OHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=155 label="OHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=155 label="OHH">, <Entry index=269 label="Cbj">]
[<Entry index=155 label="OHH">, <Entry index=168 label="OjCs">]
[<Entry index=155 label="OHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=155 label="OHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=155 label="OHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=155 label="OHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=155 label="OHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=155 label="OHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=155 label="OHH">, <Entry index=175 label="Cj">]
[<Entry index=155 label="OHH">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=155 label="OHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=155 label="OHH">, <Entry index=163 label="Hrad">]
[<Entry index=155 label="OHH">, <Entry index=335 label="C_quartetR">]
[<Entry index=155 label="OHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=155 label="OHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=155 label="OHH">, <Entry index=172 label="OjO">]
[<Entry index=155 label="OHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=177 label="Cs_methyl">]
""",
)

entry(
    index = 157,
    label = "OCH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 O u0 {2,S} {3,S}
2 *2 H u0 {1,S}
3    C ux {1,S}
""",
    kinetics = DistanceData(distances={}),
    shortDesc = u"""Fitted to 119 distances.
""",
)

entry(
    index = 158,
    label = "O/Cs/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 O  u0 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    Cs ux {1,S}
""",
    kinetics = DistanceData(distances={}),
    shortDesc = u"""Fitted to 113 distances.
""",
)

entry(
    index = 159,
    label = "O/Cd/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 O  u0 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    Cd ux {1,S}
""",
    kinetics = DistanceData(distances={}),
    shortDesc = u"""Fitted to 3 distances.
""",
)

entry(
    index = 160,
    label = "O/Ct/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 O  u0 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 161,
    label = "O/Cb/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 O  u0 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
    shortDesc = u"""Fitted to 1 distances.
""",
)

entry(
    index = 162,
    label = "OOH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *1 O u0 {2,S} {3,S}
2 *2 H u0 {1,S}
3    O ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.143199, 'd13': -0.077678, 'd23': 0.066957},
        uncertainties = {'d12': 0.126664, 'd13': 0.111413, 'd23': 0.207908},
    ),
    shortDesc = u"""Fitted to 114 distances.
""",
    longDesc = 
u"""
[<Entry index=161 label="OOH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=161 label="OOH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=161 label="OOH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=161 label="OOH">, <Entry index=184 label="CsjOH2">]
[<Entry index=161 label="OOH">, <Entry index=186 label="CsjCCH">]
[<Entry index=161 label="OOH">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=161 label="OOH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=161 label="OOH">, <Entry index=169 label="OjCd">]
[<Entry index=161 label="OOH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=161 label="OOH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=161 label="OOH">, <Entry index=175 label="Cj">]
[<Entry index=161 label="OOH">, <Entry index=172 label="OjO">]
[<Entry index=161 label="OOH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=161 label="OOH">, <Entry index=163 label="Hrad">]
[<Entry index=161 label="OOH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=161 label="OOH">, <Entry index=179 label="CsjCH2">]
[<Entry index=161 label="OOH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=161 label="OOH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=161 label="OOH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=161 label="OOH">, <Entry index=168 label="OjCs">]
""",
)

entry(
    index = 163,
    label = "Hrad",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 H u1
""",
    kinetics = DistanceData(
        distances = {'d12': -0.0269, 'd13': -0.33978, 'd23': -0.323775},
        uncertainties = {'d12': 0.155781, 'd13': 0.142754, 'd23': 0.126757},
    ),
    shortDesc = u"""Fitted to 105 distances.
""",
    longDesc = 
u"""
[<Entry index=5 label="C_methane">, <Entry index=163 label="Hrad">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=163 label="Hrad">]
[<Entry index=14 label="CsCCHH">, <Entry index=163 label="Hrad">]
[<Entry index=12 label="CsOHHH">, <Entry index=163 label="Hrad">]
[<Entry index=7 label="CsCHHH">, <Entry index=163 label="Hrad">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=163 label="Hrad">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=163 label="Hrad">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=163 label="Hrad">]
[<Entry index=2 label="C_H">, <Entry index=163 label="Hrad">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=163 label="Hrad">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=163 label="Hrad">]
[<Entry index=155 label="OHH">, <Entry index=163 label="Hrad">]
[<Entry index=71 label="C_methyl">, <Entry index=163 label="Hrad">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=151 label="Cb_H">, <Entry index=163 label="Hrad">]
[<Entry index=161 label="OOH">, <Entry index=163 label="Hrad">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=163 label="Hrad">]
[<Entry index=153 label="OradH">, <Entry index=163 label="Hrad">]
[<Entry index=25 label="CsCOHH">, <Entry index=163 label="Hrad">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=163 label="Hrad">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=163 label="Hrad">]
""",
)

entry(
    index = 164,
    label = "Orad",
    group = "OR{OjR, O_atom_triplet}",
    kinetics = DistanceData(
        distances = {'d12': -0.027946, 'd13': -0.090315, 'd23': -0.058485},
        uncertainties = {'d12': 0.177189, 'd13': 0.167298, 'd23': 0.166201},
    ),
    shortDesc = u"""Fitted to 651 distances.
""",
    longDesc = 
u"""
[<Entry index=153 label="OradH">, <Entry index=168 label="OjCs">]
[<Entry index=15 label="C/H2/Cs/Cs">, <Entry index=172 label="OjO">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=166 label="OjH">]
[<Entry index=1 label="H2">, <Entry index=168 label="OjCs">]
[<Entry index=12 label="CsOHHH">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=172 label="OjO">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=166 label="OjH">]
[<Entry index=5 label="C_methane">, <Entry index=166 label="OjH">]
[<Entry index=12 label="CsOHHH">, <Entry index=167 label="OjC">]
[<Entry index=7 label="CsCHHH">, <Entry index=166 label="OjH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=151 label="Cb_H">, <Entry index=172 label="OjO">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=172 label="OjO">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=166 label="OjH">]
[<Entry index=25 label="CsCOHH">, <Entry index=172 label="OjO">]
[<Entry index=5 label="C_methane">, <Entry index=168 label="OjCs">]
[<Entry index=155 label="OHH">, <Entry index=168 label="OjCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=14 label="CsCCHH">, <Entry index=166 label="OjH">]
[<Entry index=2 label="C_H">, <Entry index=168 label="OjCs">]
[<Entry index=14 label="CsCCHH">, <Entry index=168 label="OjCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=166 label="OjH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=168 label="OjCs">]
[<Entry index=151 label="Cb_H">, <Entry index=166 label="OjH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=172 label="OjO">]
[<Entry index=2 label="C_H">, <Entry index=167 label="OjC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=168 label="OjCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=172 label="OjO">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=71 label="C_methyl">, <Entry index=172 label="OjO">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=172 label="OjO">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=172 label="OjO">]
[<Entry index=4 label="Csnorad_H">, <Entry index=168 label="OjCs">]
[<Entry index=5 label="C_methane">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=1 label="H2">, <Entry index=172 label="OjO">]
[<Entry index=7 label="CsCHHH">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=2 label="C_H">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=71 label="C_methyl">, <Entry index=166 label="OjH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=168 label="OjCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=166 label="OjH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=166 label="OjH">]
[<Entry index=14 label="CsCCHH">, <Entry index=172 label="OjO">]
[<Entry index=161 label="OOH">, <Entry index=169 label="OjCd">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=172 label="OjO">]
[<Entry index=12 label="CsOHHH">, <Entry index=168 label="OjCs">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=172 label="OjO">]
[<Entry index=153 label="OradH">, <Entry index=172 label="OjO">]
[<Entry index=12 label="CsOHHH">, <Entry index=172 label="OjO">]
[<Entry index=161 label="OOH">, <Entry index=172 label="OjO">]
[<Entry index=1 label="H2">, <Entry index=169 label="OjCd">]
[<Entry index=7 label="CsCHHH">, <Entry index=168 label="OjCs">]
[<Entry index=71 label="C_methyl">, <Entry index=168 label="OjCs">]
[<Entry index=154 label="ORH">, <Entry index=172 label="OjO">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=155 label="OHH">, <Entry index=172 label="OjO">]
[<Entry index=1 label="H2">, <Entry index=166 label="OjH">]
[<Entry index=5 label="C_methane">, <Entry index=169 label="OjCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=166 label="OjH">]
[<Entry index=7 label="CsCHHH">, <Entry index=172 label="OjO">]
[<Entry index=1 label="H2">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=166 label="OjH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=172 label="OjO">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=168 label="OjCs">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=166 label="OjH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=2 label="C_H">, <Entry index=172 label="OjO">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=166 label="OjH">]
[<Entry index=2 label="C_H">, <Entry index=166 label="OjH">]
[<Entry index=161 label="OOH">, <Entry index=168 label="OjCs">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=172 label="OjO">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=166 label="OjH">]
[<Entry index=5 label="C_methane">, <Entry index=172 label="OjO">]
""",
)

entry(
    index = 165,
    label = "OjR",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 O u1 {2,S}
2    R ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.024349, 'd13': -0.089706, 'd23': -0.061012},
        uncertainties = {'d12': 0.18073, 'd13': 0.166315, 'd23': 0.164896},
    ),
    shortDesc = u"""Fitted to 609 distances.
""",
    longDesc = 
u"""
[<Entry index=153 label="OradH">, <Entry index=168 label="OjCs">]
[<Entry index=15 label="C/H2/Cs/Cs">, <Entry index=172 label="OjO">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=166 label="OjH">]
[<Entry index=1 label="H2">, <Entry index=168 label="OjCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=172 label="OjO">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=166 label="OjH">]
[<Entry index=5 label="C_methane">, <Entry index=166 label="OjH">]
[<Entry index=12 label="CsOHHH">, <Entry index=167 label="OjC">]
[<Entry index=7 label="CsCHHH">, <Entry index=166 label="OjH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=151 label="Cb_H">, <Entry index=172 label="OjO">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=172 label="OjO">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=166 label="OjH">]
[<Entry index=25 label="CsCOHH">, <Entry index=172 label="OjO">]
[<Entry index=5 label="C_methane">, <Entry index=168 label="OjCs">]
[<Entry index=155 label="OHH">, <Entry index=168 label="OjCs">]
[<Entry index=14 label="CsCCHH">, <Entry index=166 label="OjH">]
[<Entry index=2 label="C_H">, <Entry index=168 label="OjCs">]
[<Entry index=155 label="OHH">, <Entry index=172 label="OjO">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=166 label="OjH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=168 label="OjCs">]
[<Entry index=151 label="Cb_H">, <Entry index=166 label="OjH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=172 label="OjO">]
[<Entry index=2 label="C_H">, <Entry index=167 label="OjC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=168 label="OjCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=172 label="OjO">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=172 label="OjO">]
[<Entry index=71 label="C_methyl">, <Entry index=172 label="OjO">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=172 label="OjO">]
[<Entry index=4 label="Csnorad_H">, <Entry index=168 label="OjCs">]
[<Entry index=12 label="CsOHHH">, <Entry index=166 label="OjH">]
[<Entry index=1 label="H2">, <Entry index=172 label="OjO">]
[<Entry index=71 label="C_methyl">, <Entry index=166 label="OjH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=168 label="OjCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=166 label="OjH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=166 label="OjH">]
[<Entry index=14 label="CsCCHH">, <Entry index=172 label="OjO">]
[<Entry index=161 label="OOH">, <Entry index=169 label="OjCd">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=172 label="OjO">]
[<Entry index=12 label="CsOHHH">, <Entry index=168 label="OjCs">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=172 label="OjO">]
[<Entry index=153 label="OradH">, <Entry index=172 label="OjO">]
[<Entry index=12 label="CsOHHH">, <Entry index=172 label="OjO">]
[<Entry index=161 label="OOH">, <Entry index=172 label="OjO">]
[<Entry index=1 label="H2">, <Entry index=169 label="OjCd">]
[<Entry index=7 label="CsCHHH">, <Entry index=168 label="OjCs">]
[<Entry index=71 label="C_methyl">, <Entry index=168 label="OjCs">]
[<Entry index=154 label="ORH">, <Entry index=172 label="OjO">]
[<Entry index=1 label="H2">, <Entry index=166 label="OjH">]
[<Entry index=5 label="C_methane">, <Entry index=169 label="OjCd">]
[<Entry index=7 label="CsCHHH">, <Entry index=172 label="OjO">]
[<Entry index=14 label="CsCCHH">, <Entry index=168 label="OjCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=166 label="OjH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=172 label="OjO">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=168 label="OjCs">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=166 label="OjH">]
[<Entry index=2 label="C_H">, <Entry index=172 label="OjO">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=166 label="OjH">]
[<Entry index=2 label="C_H">, <Entry index=166 label="OjH">]
[<Entry index=161 label="OOH">, <Entry index=168 label="OjCs">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=172 label="OjO">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=166 label="OjH">]
[<Entry index=5 label="C_methane">, <Entry index=172 label="OjO">]
""",
)

entry(
    index = 166,
    label = "OjH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 O u1 {2,S}
2    H u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.15436, 'd13': -0.040208, 'd23': 0.137861},
        uncertainties = {'d12': 0.080459, 'd13': 0.217096, 'd23': 0.190701},
    ),
    shortDesc = u"""Fitted to 74 distances.
""",
    longDesc = 
u"""
[<Entry index=71 label="C_methyl">, <Entry index=166 label="OjH">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=166 label="OjH">]
[<Entry index=7 label="CsCHHH">, <Entry index=166 label="OjH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=166 label="OjH">]
[<Entry index=151 label="Cb_H">, <Entry index=166 label="OjH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=166 label="OjH">]
[<Entry index=5 label="C_methane">, <Entry index=166 label="OjH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=166 label="OjH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=166 label="OjH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=166 label="OjH">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=166 label="OjH">]
[<Entry index=14 label="CsCCHH">, <Entry index=166 label="OjH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=166 label="OjH">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=166 label="OjH">]
[<Entry index=1 label="H2">, <Entry index=166 label="OjH">]
[<Entry index=2 label="C_H">, <Entry index=166 label="OjH">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=166 label="OjH">]
[<Entry index=12 label="CsOHHH">, <Entry index=166 label="OjH">]
""",
)

entry(
    index = 167,
    label = "OjC",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 O u1 {2,S}
2    C ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.109348, 'd13': -0.116462, 'd23': -0.002109},
        uncertainties = {'d12': 0.110782, 'd13': 0.205949, 'd23': 0.189623},
    ),
    shortDesc = u"""Fitted to 118 distances.
""",
    longDesc = 
u"""
[<Entry index=1 label="H2">, <Entry index=168 label="OjCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=168 label="OjCs">]
[<Entry index=12 label="CsOHHH">, <Entry index=167 label="OjC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=5 label="C_methane">, <Entry index=168 label="OjCs">]
[<Entry index=155 label="OHH">, <Entry index=168 label="OjCs">]
[<Entry index=2 label="C_H">, <Entry index=168 label="OjCs">]
[<Entry index=5 label="C_methane">, <Entry index=169 label="OjCd">]
[<Entry index=2 label="C_H">, <Entry index=167 label="OjC">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=4 label="Csnorad_H">, <Entry index=168 label="OjCs">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=168 label="OjCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=168 label="OjCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=168 label="OjCs">]
[<Entry index=161 label="OOH">, <Entry index=169 label="OjCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=168 label="OjCs">]
[<Entry index=153 label="OradH">, <Entry index=168 label="OjCs">]
[<Entry index=1 label="H2">, <Entry index=169 label="OjCd">]
[<Entry index=7 label="CsCHHH">, <Entry index=168 label="OjCs">]
[<Entry index=71 label="C_methyl">, <Entry index=168 label="OjCs">]
[<Entry index=14 label="CsCCHH">, <Entry index=168 label="OjCs">]
[<Entry index=161 label="OOH">, <Entry index=168 label="OjCs">]
""",
)

entry(
    index = 168,
    label = "OjCs",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 O  u1 {2,S}
2    Cs ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.108557, 'd13': -0.12086, 'd23': -0.007504},
        uncertainties = {'d12': 0.10667, 'd13': 0.18317, 'd23': 0.15138},
    ),
    shortDesc = u"""Fitted to 112 distances.
""",
    longDesc = 
u"""
[<Entry index=71 label="C_methyl">, <Entry index=168 label="OjCs">]
[<Entry index=1 label="H2">, <Entry index=168 label="OjCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=168 label="OjCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=168 label="OjCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=168 label="OjCs">]
[<Entry index=14 label="CsCCHH">, <Entry index=168 label="OjCs">]
[<Entry index=2 label="C_H">, <Entry index=168 label="OjCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=5 label="C_methane">, <Entry index=168 label="OjCs">]
[<Entry index=12 label="CsOHHH">, <Entry index=168 label="OjCs">]
[<Entry index=155 label="OHH">, <Entry index=168 label="OjCs">]
[<Entry index=153 label="OradH">, <Entry index=168 label="OjCs">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=168 label="OjCs">]
[<Entry index=4 label="Csnorad_H">, <Entry index=168 label="OjCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=168 label="OjCs">]
[<Entry index=7 label="CsCHHH">, <Entry index=168 label="OjCs">]
[<Entry index=161 label="OOH">, <Entry index=168 label="OjCs">]
""",
)

entry(
    index = 169,
    label = "OjCd",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 O  u1 {2,S}
2    Cd ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.075299, 'd13': -0.147585, 'd23': -0.074439},
        uncertainties = {'d12': 0.240134, 'd13': 0.137748, 'd23': 0.117701},
    ),
    shortDesc = u"""Fitted to 3 distances.
""",
    longDesc = 
u"""
[<Entry index=1 label="H2">, <Entry index=169 label="OjCd">]
[<Entry index=5 label="C_methane">, <Entry index=169 label="OjCd">]
[<Entry index=161 label="OOH">, <Entry index=169 label="OjCd">]
""",
)

entry(
    index = 170,
    label = "OjCt",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 O  u1 {2,S}
2    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 171,
    label = "OjCb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 O  u1 {2,S}
2    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
    shortDesc = u"""Fitted to 1 distances.
""",
)

entry(
    index = 172,
    label = "OjO",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 O u1 {2,S}
2    O ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.019609, 'd13': -0.092246, 'd23': -0.111175},
        uncertainties = {'d12': 0.208138, 'd13': 0.143889, 'd23': 0.153763},
    ),
    shortDesc = u"""Fitted to 417 distances.
""",
    longDesc = 
u"""
[<Entry index=9 label="C/H3/Cd">, <Entry index=172 label="OjO">]
[<Entry index=151 label="Cb_H">, <Entry index=172 label="OjO">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=172 label="OjO">]
[<Entry index=25 label="CsCOHH">, <Entry index=172 label="OjO">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=172 label="OjO">]
[<Entry index=15 label="C/H2/Cs/Cs">, <Entry index=172 label="OjO">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=172 label="OjO">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=172 label="OjO">]
[<Entry index=71 label="C_methyl">, <Entry index=172 label="OjO">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=172 label="OjO">]
[<Entry index=1 label="H2">, <Entry index=172 label="OjO">]
[<Entry index=14 label="CsCCHH">, <Entry index=172 label="OjO">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=172 label="OjO">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=172 label="OjO">]
[<Entry index=153 label="OradH">, <Entry index=172 label="OjO">]
[<Entry index=12 label="CsOHHH">, <Entry index=172 label="OjO">]
[<Entry index=161 label="OOH">, <Entry index=172 label="OjO">]
[<Entry index=154 label="ORH">, <Entry index=172 label="OjO">]
[<Entry index=7 label="CsCHHH">, <Entry index=172 label="OjO">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=172 label="OjO">]
[<Entry index=155 label="OHH">, <Entry index=172 label="OjO">]
[<Entry index=2 label="C_H">, <Entry index=172 label="OjO">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=172 label="OjO">]
[<Entry index=5 label="C_methane">, <Entry index=172 label="OjO">]
""",
)

entry(
    index = 173,
    label = "O_atom_triplet",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 O u2
""",
    kinetics = DistanceData(
        distances = {'d12': -0.078337, 'd13': -0.098843, 'd23': -0.023085},
        uncertainties = {'d12': 0.119453, 'd13': 0.188641, 'd23': 0.191885},
    ),
    shortDesc = u"""Fitted to 42 distances.
""",
    longDesc = 
u"""
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=7 label="CsCHHH">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=12 label="CsOHHH">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=1 label="H2">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=2 label="C_H">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=173 label="O_atom_triplet">]
[<Entry index=5 label="C_methane">, <Entry index=173 label="O_atom_triplet">]
""",
)

entry(
    index = 174,
    label = "Crad",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 C u[1,2,3,4]
""",
    kinetics = DistanceData(
        distances = {'d12': 0.013496, 'd13': 0.059843, 'd23': 0.045567},
        uncertainties = {'d12': 0.15759, 'd13': 0.13968, 'd23': 0.141434},
    ),
    shortDesc = u"""Fitted to 1735 distances.
""",
    longDesc = 
u"""
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=1 label="H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=151 label="Cb_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=268 label="CtjC">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=150 label="Ct_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=155 label="OHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=269 label="Cbj">]
[<Entry index=5 label="C_methane">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=184 label="CsjOH2">]
[<Entry index=4 label="Csnorad_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=2 label="C_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=155 label="OHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=161 label="OOH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=269 label="Cbj">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=204 label="CsjCCC">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=151 label="Cb_H">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=155 label="OHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=5 label="C_methane">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=2 label="C_H">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=155 label="OHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=5 label="C_methane">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=161 label="OOH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=155 label="OHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=1 label="H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=71 label="C_methyl">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=5 label="C_methane">, <Entry index=175 label="Cj">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=161 label="OOH">, <Entry index=186 label="CsjCCH">]
[<Entry index=151 label="Cb_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=71 label="C_methyl">, <Entry index=204 label="CsjCCC">]
[<Entry index=12 label="CsOHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=204 label="CsjCCC">]
[<Entry index=12 label="CsOHHH">, <Entry index=269 label="Cbj">]
[<Entry index=71 label="C_methyl">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=153 label="OradH">, <Entry index=186 label="CsjCCH">]
[<Entry index=1 label="H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=161 label="OOH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=1 label="H2">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=5 label="C_methane">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=12 label="CsOHHH">, <Entry index=175 label="Cj">]
[<Entry index=161 label="OOH">, <Entry index=184 label="CsjOH2">]
[<Entry index=150 label="Ct_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=161 label="OOH">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=71 label="C_methyl">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=7 label="CsCHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=153 label="OradH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=1 label="H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=175 label="Cj">]
[<Entry index=154 label="ORH">, <Entry index=184 label="CsjOH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=153 label="OradH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=151 label="Cb_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=177 label="Cs_methyl">]
[<Entry index=155 label="OHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=150 label="Ct_H">, <Entry index=268 label="CtjC">]
[<Entry index=150 label="Ct_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=175 label="Cj">]
[<Entry index=7 label="CsCHHH">, <Entry index=175 label="Cj">]
[<Entry index=1 label="H2">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=1 label="H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=2 label="C_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=150 label="Ct_H">, <Entry index=269 label="Cbj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=1 label="H2">, <Entry index=175 label="Cj">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=150 label="Ct_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=179 label="CsjCH2">]
[<Entry index=127 label="Cd_Cds/Cd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=2 label="C_H">, <Entry index=175 label="Cj">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=204 label="CsjCCC">]
[<Entry index=150 label="Ct_H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=161 label="OOH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=174 label="Crad">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=161 label="OOH">, <Entry index=179 label="CsjCH2">]
[<Entry index=161 label="OOH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=150 label="Ct_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=153 label="OradH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=1 label="H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=161 label="OOH">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=71 label="C_methyl">, <Entry index=175 label="Cj">]
[<Entry index=155 label="OHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=153 label="OradH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=151 label="Cb_H">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=175 label="Cj">]
[<Entry index=7 label="CsCHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=175 label="Cj">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=175 label="Cj">]
[<Entry index=155 label="OHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=1 label="H2">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=12 label="CsOHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=1 label="H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=153 label="OradH">, <Entry index=204 label="CsjCCC">]
[<Entry index=151 label="Cb_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=2 label="C_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=268 label="CtjC">]
[<Entry index=121 label="Cd_C/R/H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=12 label="CsOHHH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=151 label="Cb_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=153 label="OradH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=5 label="C_methane">, <Entry index=186 label="CsjCCH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=4 label="Csnorad_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=151 label="Cb_H">, <Entry index=175 label="Cj">]
[<Entry index=161 label="OOH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=7 label="CsCHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=175 label="Cj">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=174 label="Crad">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=155 label="OHH">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=161 label="OOH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=2 label="C_H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=4 label="Csnorad_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=174 label="Crad">]
[<Entry index=2 label="C_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=153 label="OradH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=2 label="C_H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=1 label="H2">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=153 label="OradH">, <Entry index=175 label="Cj">]
[<Entry index=153 label="OradH">, <Entry index=179 label="CsjCH2">]
[<Entry index=151 label="Cb_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=153 label="OradH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=12 label="CsOHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=4 label="Csnorad_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=151 label="Cb_H">, <Entry index=269 label="Cbj">]
[<Entry index=155 label="OHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=204 label="CsjCCC">]
[<Entry index=5 label="C_methane">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=186 label="CsjCCH">]
[<Entry index=1 label="H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=153 label="OradH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=1 label="H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=5 label="C_methane">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=151 label="Cb_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=151 label="Cb_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=153 label="OradH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=1 label="H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=4 label="Csnorad_H">, <Entry index=175 label="Cj">]
[<Entry index=12 label="CsOHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=180 label="Csj/Cs/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=179 label="CsjCH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=184 label="CsjOH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=268 label="CtjC">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=1 label="H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=161 label="OOH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=1 label="H2">, <Entry index=204 label="CsjCCC">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=5 label="C_methane">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=12 label="CsOHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=186 label="CsjCCH">]
[<Entry index=5 label="C_methane">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=268 label="CtjC">]
[<Entry index=155 label="OHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=155 label="OHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=2 label="C_H">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=269 label="Cbj">]
[<Entry index=2 label="C_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=1 label="H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=5 label="C_methane">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=155 label="OHH">, <Entry index=335 label="C_quartetR">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=161 label="OOH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=71 label="C_methyl">, <Entry index=186 label="CsjCCH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=161 label="OOH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=1 label="H2">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=5 label="C_methane">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=71 label="C_methyl">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=155 label="OHH">, <Entry index=269 label="Cbj">]
[<Entry index=150 label="Ct_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=7 label="CsCHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=1 label="H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=2 label="C_H">, <Entry index=269 label="Cbj">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=161 label="OOH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=2 label="C_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=14 label="CsCCHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=269 label="Cbj">]
[<Entry index=150 label="Ct_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=161 label="OOH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=175 label="Cj">]
[<Entry index=4 label="Csnorad_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=151 label="Cb_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=1 label="H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=4 label="Csnorad_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=155 label="OHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=2 label="C_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=150 label="Ct_H">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=269 label="Cbj">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=1 label="H2">, <Entry index=186 label="CsjCCH">]
[<Entry index=1 label="H2">, <Entry index=335 label="C_quartetR">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=5 label="C_methane">, <Entry index=184 label="CsjOH2">]
[<Entry index=153 label="OradH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=2 label="C_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=5 label="C_methane">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=151 label="Cb_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=12 label="CsOHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=5 label="C_methane">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=186 label="CsjCCH">]
[<Entry index=7 label="CsCHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=177 label="Cs_methyl">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=7 label="CsCHHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=150 label="Ct_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=153 label="OradH">, <Entry index=184 label="CsjOH2">]
[<Entry index=155 label="OHH">, <Entry index=175 label="Cj">]
[<Entry index=5 label="C_methane">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=2 label="C_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=2 label="C_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=1 label="H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=25 label="CsCOHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=12 label="CsOHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=2 label="C_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=335 label="C_quartetR">]
[<Entry index=7 label="CsCHHH">, <Entry index=269 label="Cbj">]
""",
)

entry(
    index = 175,
    label = "Cj",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 C u1
""",
    kinetics = DistanceData(
        distances = {'d12': 0.013477, 'd13': 0.060538, 'd23': 0.045679},
        uncertainties = {'d12': 0.157238, 'd13': 0.139704, 'd23': 0.14024},
    ),
    shortDesc = u"""Fitted to 1668 distances.
""",
    longDesc = 
u"""
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=161 label="OOH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=1 label="H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=1 label="H2">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=4 label="Csnorad_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=5 label="C_methane">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=184 label="CsjOH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=71 label="C_methyl">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=155 label="OHH">, <Entry index=269 label="Cbj">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=1 label="H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=151 label="Cb_H">, <Entry index=175 label="Cj">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=151 label="Cb_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=150 label="Ct_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=175 label="Cj">]
[<Entry index=127 label="Cd_Cds/Cd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=161 label="OOH">, <Entry index=179 label="CsjCH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=151 label="Cb_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=268 label="CtjC">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=150 label="Ct_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=175 label="Cj">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=2 label="C_H">, <Entry index=269 label="Cbj">]
[<Entry index=155 label="OHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=1 label="H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=269 label="Cbj">]
[<Entry index=4 label="Csnorad_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=153 label="OradH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=153 label="OradH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=161 label="OOH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=2 label="C_H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=2 label="C_H">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=150 label="Ct_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=4 label="Csnorad_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=2 label="C_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=2 label="C_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=2 label="C_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=71 label="C_methyl">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=155 label="OHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=2 label="C_H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=5 label="C_methane">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=161 label="OOH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=269 label="Cbj">]
[<Entry index=5 label="C_methane">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=269 label="Cbj">]
[<Entry index=150 label="Ct_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=161 label="OOH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=153 label="OradH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=1 label="H2">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=12 label="CsOHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=4 label="Csnorad_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=1 label="H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=151 label="Cb_H">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=153 label="OradH">, <Entry index=175 label="Cj">]
[<Entry index=1 label="H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=153 label="OradH">, <Entry index=179 label="CsjCH2">]
[<Entry index=2 label="C_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=177 label="Cs_methyl">]
[<Entry index=151 label="Cb_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=12 label="CsOHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=5 label="C_methane">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=4 label="Csnorad_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=5 label="C_methane">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=151 label="Cb_H">, <Entry index=269 label="Cbj">]
[<Entry index=155 label="OHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=5 label="C_methane">, <Entry index=186 label="CsjCCH">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=204 label="CsjCCC">]
[<Entry index=14 label="CsCCHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=155 label="OHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=204 label="CsjCCC">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=2 label="C_H">, <Entry index=175 label="Cj">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=161 label="OOH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=1 label="H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=153 label="OradH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=5 label="C_methane">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=151 label="Cb_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=1 label="H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=5 label="C_methane">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=4 label="Csnorad_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=7 label="CsCHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=150 label="Ct_H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=186 label="CsjCCH">]
[<Entry index=71 label="C_methyl">, <Entry index=175 label="Cj">]
[<Entry index=155 label="OHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=161 label="OOH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=153 label="OradH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=180 label="Csj/Cs/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=5 label="C_methane">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=2 label="C_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=269 label="Cbj">]
[<Entry index=1 label="H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=161 label="OOH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=1 label="H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=161 label="OOH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=151 label="Cb_H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=150 label="Ct_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=5 label="C_methane">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=184 label="CsjOH2">]
[<Entry index=155 label="OHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=151 label="Cb_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=7 label="CsCHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=204 label="CsjCCC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=1 label="H2">, <Entry index=186 label="CsjCCH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=5 label="C_methane">, <Entry index=175 label="Cj">]
[<Entry index=2 label="C_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=161 label="OOH">, <Entry index=186 label="CsjCCH">]
[<Entry index=153 label="OradH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=1 label="H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=151 label="Cb_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=5 label="C_methane">, <Entry index=179 label="CsjCH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=5 label="C_methane">, <Entry index=184 label="CsjOH2">]
[<Entry index=71 label="C_methyl">, <Entry index=204 label="CsjCCC">]
[<Entry index=4 label="Csnorad_H">, <Entry index=175 label="Cj">]
[<Entry index=12 label="CsOHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=177 label="Cs_methyl">]
[<Entry index=153 label="OradH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=5 label="C_methane">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=186 label="CsjCCH">]
[<Entry index=12 label="CsOHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=179 label="CsjCH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=184 label="CsjOH2">]
[<Entry index=5 label="C_methane">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=161 label="OOH">, <Entry index=175 label="Cj">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=268 label="CtjC">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=204 label="CsjCCC">]
[<Entry index=12 label="CsOHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=12 label="CsOHHH">, <Entry index=269 label="Cbj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=150 label="Ct_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=153 label="OradH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=155 label="OHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=161 label="OOH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=1 label="H2">, <Entry index=204 label="CsjCCC">]
[<Entry index=5 label="C_methane">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=155 label="OHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=175 label="Cj">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=175 label="Cj">]
[<Entry index=155 label="OHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=7 label="CsCHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=186 label="CsjCCH">]
[<Entry index=12 label="CsOHHH">, <Entry index=175 label="Cj">]
[<Entry index=5 label="C_methane">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=1 label="H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=186 label="CsjCCH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=268 label="CtjC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=155 label="OHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=175 label="Cj">]
[<Entry index=155 label="OHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=175 label="Cj">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=175 label="Cj">]
[<Entry index=71 label="C_methyl">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=155 label="OHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=2 label="C_H">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=150 label="Ct_H">, <Entry index=269 label="Cbj">]
[<Entry index=12 label="CsOHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=269 label="Cbj">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=150 label="Ct_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=1 label="H2">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=2 label="C_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=7 label="CsCHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=153 label="OradH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=1 label="H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=175 label="Cj">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=175 label="Cj">]
[<Entry index=5 label="C_methane">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=1 label="H2">, <Entry index=175 label="Cj">]
[<Entry index=153 label="OradH">, <Entry index=204 label="CsjCCC">]
[<Entry index=151 label="Cb_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=154 label="ORH">, <Entry index=184 label="CsjOH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=1 label="H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=2 label="C_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=153 label="OradH">, <Entry index=186 label="CsjCCH">]
[<Entry index=151 label="Cb_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=151 label="Cb_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=153 label="OradH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=175 label="Cj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=268 label="CtjC">]
[<Entry index=121 label="Cd_C/R/H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=2 label="C_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=1 label="H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=151 label="Cb_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=2 label="C_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=268 label="CtjC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=161 label="OOH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=2 label="C_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=161 label="OOH">, <Entry index=184 label="CsjOH2">]
[<Entry index=25 label="CsCOHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=269 label="Cbj">]
[<Entry index=153 label="OradH">, <Entry index=191 label="Csj/Cd/Cd/H">]
""",
)

entry(
    index = 176,
    label = "Csj",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1
""",
    kinetics = DistanceData(
        distances = {'d12': 0.025684, 'd13': 0.057431, 'd23': 0.028584},
        uncertainties = {'d12': 0.145986, 'd13': 0.131434, 'd23': 0.129109},
    ),
    shortDesc = u"""Fitted to 1112 distances.
""",
    longDesc = 
u"""
[<Entry index=12 label="CsOHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=161 label="OOH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=1 label="H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=179 label="CsjCH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=4 label="Csnorad_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=71 label="C_methyl">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=1 label="H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=150 label="Ct_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=204 label="CsjCCC">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=151 label="Cb_H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=151 label="Cb_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=7 label="CsCHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=150 label="Ct_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=155 label="OHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=161 label="OOH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=2 label="C_H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=2 label="C_H">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=4 label="Csnorad_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=2 label="C_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=2 label="C_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=2 label="C_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=14 label="CsCCHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=2 label="C_H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=5 label="C_methane">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=153 label="OradH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=1 label="H2">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=4 label="Csnorad_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=153 label="OradH">, <Entry index=184 label="CsjOH2">]
[<Entry index=1 label="H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=153 label="OradH">, <Entry index=179 label="CsjCH2">]
[<Entry index=2 label="C_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=151 label="Cb_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=127 label="Cd_Cds/Cd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=155 label="OHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=204 label="CsjCCC">]
[<Entry index=155 label="OHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=161 label="OOH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=1 label="H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=153 label="OradH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=1 label="H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=186 label="CsjCCH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=204 label="CsjCCC">]
[<Entry index=150 label="Ct_H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=161 label="OOH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=180 label="Csj/Cs/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=5 label="C_methane">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=1 label="H2">, <Entry index=204 label="CsjCCC">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=1 label="H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=161 label="OOH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=150 label="Ct_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=5 label="C_methane">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=2 label="C_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=151 label="Cb_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=204 label="CsjCCC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=1 label="H2">, <Entry index=186 label="CsjCCH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=2 label="C_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=161 label="OOH">, <Entry index=186 label="CsjCCH">]
[<Entry index=153 label="OradH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=1 label="H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=151 label="Cb_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=5 label="C_methane">, <Entry index=179 label="CsjCH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=5 label="C_methane">, <Entry index=184 label="CsjOH2">]
[<Entry index=71 label="C_methyl">, <Entry index=204 label="CsjCCC">]
[<Entry index=12 label="CsOHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=177 label="Cs_methyl">]
[<Entry index=153 label="OradH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=179 label="CsjCH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=184 label="CsjOH2">]
[<Entry index=5 label="C_methane">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=204 label="CsjCCC">]
[<Entry index=12 label="CsOHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=150 label="Ct_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=155 label="OHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=161 label="OOH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=1 label="H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=5 label="C_methane">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=161 label="OOH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=5 label="C_methane">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=12 label="CsOHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=186 label="CsjCCH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=186 label="CsjCCH">]
[<Entry index=1 label="H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=186 label="CsjCCH">]
[<Entry index=161 label="OOH">, <Entry index=184 label="CsjOH2">]
[<Entry index=155 label="OHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=177 label="Cs_methyl">]
[<Entry index=155 label="OHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=155 label="OHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=2 label="C_H">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=71 label="C_methyl">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=150 label="Ct_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=1 label="H2">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=2 label="C_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=153 label="OradH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=1 label="H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=5 label="C_methane">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=154 label="ORH">, <Entry index=184 label="CsjOH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=184 label="CsjOH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=186 label="CsjCCH">]
[<Entry index=151 label="Cb_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=5 label="C_methane">, <Entry index=186 label="CsjCCH">]
[<Entry index=153 label="OradH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=121 label="Cd_C/R/H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=4 label="Csnorad_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=161 label="OOH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=151 label="Cb_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=2 label="C_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=7 label="CsCHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=2 label="C_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=25 label="CsCOHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=71 label="C_methyl">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=191 label="Csj/Cd/Cd/H">]
""",
)

entry(
    index = 177,
    label = "Cs_methyl",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.04032, 'd13': 0.073629, 'd23': 0.102813},
        uncertainties = {'d12': 0.112393, 'd13': 0.145185, 'd23': 0.117782},
    ),
    shortDesc = u"""Fitted to 82 distances.
""",
    longDesc = 
u"""
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=7 label="CsCHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=151 label="Cb_H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=1 label="H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=177 label="Cs_methyl">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=26 label="C/H2/Cs/O">, <Entry index=177 label="Cs_methyl">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=177 label="Cs_methyl">]
[<Entry index=153 label="OradH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=155 label="OHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=14 label="CsCCHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=12 label="CsOHHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=177 label="Cs_methyl">]
[<Entry index=161 label="OOH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=25 label="CsCOHH">, <Entry index=177 label="Cs_methyl">]
[<Entry index=2 label="C_H">, <Entry index=177 label="Cs_methyl">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=177 label="Cs_methyl">]
""",
)

entry(
    index = 178,
    label = "CsjRH2",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs  u1 {2,S} {3,S} {4,S}
2    H   u0 {1,S}
3    H   u0 {1,S}
4    R!H ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.012695, 'd13': 0.052917, 'd23': 0.036101},
        uncertainties = {'d12': 0.146304, 'd13': 0.139006, 'd23': 0.138485},
    ),
    shortDesc = u"""Fitted to 637 distances.
""",
    longDesc = 
u"""
[<Entry index=12 label="CsOHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=1 label="H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=1 label="H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=161 label="OOH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=155 label="OHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=1 label="H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=5 label="C_methane">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=4 label="Csnorad_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=5 label="C_methane">, <Entry index=179 label="CsjCH2">]
[<Entry index=153 label="OradH">, <Entry index=179 label="CsjCH2">]
[<Entry index=2 label="C_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=2 label="C_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=155 label="OHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=1 label="H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=161 label="OOH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=180 label="Csj/Cs/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=153 label="OradH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=151 label="Cb_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=151 label="Cb_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=151 label="Cb_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=5 label="C_methane">, <Entry index=184 label="CsjOH2">]
[<Entry index=153 label="OradH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=179 label="CsjCH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=184 label="CsjOH2">]
[<Entry index=150 label="Ct_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=1 label="H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=161 label="OOH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=161 label="OOH">, <Entry index=179 label="CsjCH2">]
[<Entry index=161 label="OOH">, <Entry index=184 label="CsjOH2">]
[<Entry index=155 label="OHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=155 label="OHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=2 label="C_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=153 label="OradH">, <Entry index=184 label="CsjOH2">]
[<Entry index=5 label="C_methane">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=5 label="C_methane">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=154 label="ORH">, <Entry index=184 label="CsjOH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=184 label="CsjOH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=121 label="Cd_C/R/H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=2 label="C_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=2 label="C_H">, <Entry index=180 label="Csj/Cs/H2">]
""",
)

entry(
    index = 179,
    label = "CsjCH2",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    H  u0 {1,S}
4    C  ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.009663, 'd13': 0.053981, 'd23': 0.038571},
        uncertainties = {'d12': 0.142898, 'd13': 0.134135, 'd23': 0.136394},
    ),
    shortDesc = u"""Fitted to 550 distances.
""",
    longDesc = 
u"""
[<Entry index=151 label="Cb_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=5 label="C_methane">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=150 label="Ct_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=4 label="Csnorad_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=5 label="C_methane">, <Entry index=179 label="CsjCH2">]
[<Entry index=1 label="H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=161 label="OOH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=153 label="OradH">, <Entry index=179 label="CsjCH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=2 label="C_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=180 label="Csj/Cs/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=155 label="OHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=155 label="OHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=5 label="C_methane">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=161 label="OOH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=1 label="H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=155 label="OHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=179 label="CsjCH2">]
[<Entry index=1 label="H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=1 label="H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=153 label="OradH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=179 label="CsjCH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=179 label="CsjCH2">]
[<Entry index=5 label="C_methane">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=179 label="CsjCH2">]
[<Entry index=161 label="OOH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=121 label="Cd_C/R/H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=161 label="OOH">, <Entry index=179 label="CsjCH2">]
[<Entry index=2 label="C_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=2 label="C_H">, <Entry index=179 label="CsjCH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=179 label="CsjCH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=153 label="OradH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=2 label="C_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=181 label="Csj/Cd/H2">]
""",
)

entry(
    index = 180,
    label = "Csj/Cs/H2",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    H  u0 {1,S}
4    Cs ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.020367, 'd13': 0.043785, 'd23': 0.056832},
        uncertainties = {'d12': 0.127062, 'd13': 0.14254, 'd23': 0.147342},
    ),
    shortDesc = u"""Fitted to 305 distances.
""",
    longDesc = 
u"""
[<Entry index=8 label="C/H3/Cs">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=4 label="Csnorad_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=161 label="OOH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=1 label="H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=155 label="OHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=108 label="Cs_tripletHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=5 label="C_methane">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=180 label="Csj/Cs/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=153 label="OradH">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=180 label="Csj/Cs/H2">]
[<Entry index=2 label="C_H">, <Entry index=180 label="Csj/Cs/H2">]
""",
)

entry(
    index = 181,
    label = "Csj/Cd/H2",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    H  u0 {1,S}
4    Cd ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.053703, 'd13': 0.0634, 'd23': 0.005333},
        uncertainties = {'d12': 0.164461, 'd13': 0.131746, 'd23': 0.114067},
    ),
    shortDesc = u"""Fitted to 156 distances.
""",
    longDesc = 
u"""
[<Entry index=121 label="Cd_C/R/H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=5 label="C_methane">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=153 label="OradH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=27 label="C/H2/Cd/O">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=155 label="OHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=2 label="C_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=1 label="H2">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=12 label="CsOHHH">, <Entry index=181 label="Csj/Cd/H2">]
[<Entry index=161 label="OOH">, <Entry index=181 label="Csj/Cd/H2">]
""",
)

entry(
    index = 182,
    label = "Csj/Ct/H2",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    H  u0 {1,S}
4    Ct u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.040376, 'd13': 0.081579, 'd23': 0.041537},
        uncertainties = {'d12': 0.216834, 'd13': 0.134118, 'd23': 0.135313},
    ),
    shortDesc = u"""Fitted to 22 distances.
""",
    longDesc = 
u"""
[<Entry index=2 label="C_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=1 label="H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=5 label="C_methane">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=151 label="Cb_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=150 label="Ct_H">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=71 label="C_methyl">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=7 label="CsCHHH">, <Entry index=182 label="Csj/Ct/H2">]
[<Entry index=161 label="OOH">, <Entry index=182 label="Csj/Ct/H2">]
""",
)

entry(
    index = 183,
    label = "Csj/Cb/H2",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    H  u0 {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 184,
    label = "CsjOH2",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    H  u0 {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.032501, 'd13': 0.045968, 'd23': 0.019964},
        uncertainties = {'d12': 0.169458, 'd13': 0.169729, 'd23': 0.154002},
    ),
    shortDesc = u"""Fitted to 87 distances.
""",
    longDesc = 
u"""
[<Entry index=2 label="C_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=153 label="OradH">, <Entry index=184 label="CsjOH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=184 label="CsjOH2">]
[<Entry index=154 label="ORH">, <Entry index=184 label="CsjOH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=184 label="CsjOH2">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=155 label="OHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=184 label="CsjOH2">]
[<Entry index=150 label="Ct_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=5 label="C_methane">, <Entry index=184 label="CsjOH2">]
[<Entry index=1 label="H2">, <Entry index=184 label="CsjOH2">]
[<Entry index=151 label="Cb_H">, <Entry index=184 label="CsjOH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=184 label="CsjOH2">]
[<Entry index=161 label="OOH">, <Entry index=184 label="CsjOH2">]
""",
)

entry(
    index = 185,
    label = "CsjRRH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs  u1 {2,S} {3,S} {4,S}
2    H   u0 {1,S}
3    R!H ux {1,S}
4    R!H ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.066937, 'd13': 0.062701, 'd23': -0.003478},
        uncertainties = {'d12': 0.150361, 'd13': 0.110083, 'd23': 0.116215},
    ),
    shortDesc = u"""Fitted to 356 distances.
""",
    longDesc = 
u"""
[<Entry index=153 label="OradH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=161 label="OOH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=4 label="Csnorad_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=186 label="CsjCCH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=186 label="CsjCCH">]
[<Entry index=1 label="H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=150 label="Ct_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=1 label="H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=153 label="OradH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=2 label="C_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=161 label="OOH">, <Entry index=186 label="CsjCCH">]
[<Entry index=155 label="OHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=151 label="Cb_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=150 label="Ct_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=127 label="Cd_Cds/Cd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=5 label="C_methane">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=1 label="H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=150 label="Ct_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=1 label="H2">, <Entry index=186 label="CsjCCH">]
[<Entry index=7 label="CsCHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=150 label="Ct_H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=1 label="H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=186 label="CsjCCH">]
[<Entry index=12 label="CsOHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=186 label="CsjCCH">]
[<Entry index=71 label="C_methyl">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=155 label="OHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=2 label="C_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=2 label="C_H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=186 label="CsjCCH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=153 label="OradH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=4 label="Csnorad_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=2 label="C_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=71 label="C_methyl">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=71 label="C_methyl">, <Entry index=186 label="CsjCCH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=151 label="Cb_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
""",
)

entry(
    index = 186,
    label = "CsjCCH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    C  ux {1,S}
4    C  ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.071204, 'd13': 0.066312, 'd23': -0.00698},
        uncertainties = {'d12': 0.149662, 'd13': 0.109503, 'd23': 0.111823},
    ),
    shortDesc = u"""Fitted to 287 distances.
""",
    longDesc = 
u"""
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=186 label="CsjCCH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=186 label="CsjCCH">]
[<Entry index=1 label="H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=150 label="Ct_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=1 label="H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=153 label="OradH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=2 label="C_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=161 label="OOH">, <Entry index=186 label="CsjCCH">]
[<Entry index=155 label="OHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=151 label="Cb_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=150 label="Ct_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=127 label="Cd_Cds/Cd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=186 label="CsjCCH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=1 label="H2">, <Entry index=186 label="CsjCCH">]
[<Entry index=7 label="CsCHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=150 label="Ct_H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=1 label="H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=186 label="CsjCCH">]
[<Entry index=71 label="C_methyl">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=155 label="OHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=2 label="C_H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=4 label="Csnorad_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=2 label="C_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=71 label="C_methyl">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=186 label="CsjCCH">]
[<Entry index=71 label="C_methyl">, <Entry index=186 label="CsjCCH">]
[<Entry index=161 label="OOH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=186 label="CsjCCH">]
[<Entry index=155 label="OHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
""",
)

entry(
    index = 187,
    label = "Csj/Cs/Cs/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cs ux {1,S}
4    Cs ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.006032, 'd13': 0.024451, 'd23': 0.025171},
        uncertainties = {'d12': 0.135678, 'd13': 0.132833, 'd23': 0.124882},
    ),
    shortDesc = u"""Fitted to 97 distances.
""",
    longDesc = 
u"""
[<Entry index=161 label="OOH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=2 label="C_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=5 label="C_methane">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=1 label="H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=150 label="Ct_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=151 label="Cb_H">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=153 label="OradH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
[<Entry index=155 label="OHH">, <Entry index=187 label="Csj/Cs/Cs/H">]
""",
)

entry(
    index = 188,
    label = "Csj/Cs/Cd/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cs ux {1,S}
4    Cd ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.100333, 'd13': 0.079552, 'd23': -0.020441},
        uncertainties = {'d12': 0.175551, 'd13': 0.103891, 'd23': 0.114786},
    ),
    shortDesc = u"""Fitted to 109 distances.
""",
    longDesc = 
u"""
[<Entry index=7 label="CsCHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=2 label="C_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=71 label="C_methyl">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=150 label="Ct_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=1 label="H2">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=4 label="Csnorad_H">, <Entry index=188 label="Csj/Cs/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=188 label="Csj/Cs/Cd/H">]
""",
)

entry(
    index = 189,
    label = "Csj/Cs/Ct/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cs ux {1,S}
4    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 190,
    label = "Csj/Cs/Cb/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cs ux {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 191,
    label = "Csj/Cd/Cd/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cd ux {1,S}
4    Cd ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.144334, 'd13': 0.112262, 'd23': -0.033452},
        uncertainties = {'d12': 0.121907, 'd13': 0.086273, 'd23': 0.08528},
    ),
    shortDesc = u"""Fitted to 59 distances.
""",
    longDesc = 
u"""
[<Entry index=10 label="C/H3/Ct">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=1 label="H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=153 label="OradH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=150 label="Ct_H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=127 label="Cd_Cds/Cd/H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=155 label="OHH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=71 label="C_methyl">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=2 label="C_H">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=161 label="OOH">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=5 label="C_methane">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=191 label="Csj/Cd/Cd/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=191 label="Csj/Cd/Cd/H">]
""",
)

entry(
    index = 192,
    label = "Csj/Cd/Ct/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cd ux {1,S}
4    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 193,
    label = "Csj/Cd/Cb/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cd ux {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 194,
    label = "Csj/Ct/Ct/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Ct u0 {1,S}
4    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 195,
    label = "Csj/Ct/Cb/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Ct u0 {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 196,
    label = "Csj/Cb/Cb/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cb u0 {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 197,
    label = "CsjCOH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    C  ux {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.047561, 'd13': 0.046305, 'd23': 0.012424},
        uncertainties = {'d12': 0.157115, 'd13': 0.115313, 'd23': 0.136142},
    ),
    shortDesc = u"""Fitted to 69 distances.
""",
    longDesc = 
u"""
[<Entry index=153 label="OradH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=155 label="OHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=1 label="H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=4 label="Csnorad_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=5 label="C_methane">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=2 label="C_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=161 label="OOH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=150 label="Ct_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=151 label="Cb_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
""",
)

entry(
    index = 198,
    label = "Csj/Cs/O/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cs ux {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.047561, 'd13': 0.046305, 'd23': 0.012424},
        uncertainties = {'d12': 0.157115, 'd13': 0.115313, 'd23': 0.136142},
    ),
    shortDesc = u"""Fitted to 69 distances.
""",
    longDesc = 
u"""
[<Entry index=153 label="OradH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=12 label="CsOHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=155 label="OHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=1 label="H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=4 label="Csnorad_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=5 label="C_methane">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=2 label="C_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=161 label="OOH">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=150 label="Ct_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=151 label="Cb_H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=198 label="Csj/Cs/O/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=198 label="Csj/Cs/O/H">]
""",
)

entry(
    index = 199,
    label = "Csj/Cd/O/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cd ux {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
    shortDesc = u"""Fitted to 1 distances.
""",
)

entry(
    index = 200,
    label = "Csj/Ct/O/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Ct u0 {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 201,
    label = "Csj/Cb/O/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cb u0 {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 202,
    label = "CsjOOH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    O  ux {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
    shortDesc = u"""Fitted to 2 distances.
""",
)

entry(
    index = 203,
    label = "CsjRRR",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs  u1 {2,S} {3,S} {4,S}
2    R!H ux {1,S}
3    R!H ux {1,S}
4    R!H ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.026887, 'd13': 0.040537, 'd23': 0.010278},
        uncertainties = {'d12': 0.180874, 'd13': 0.170619, 'd23': 0.118981},
    ),
    shortDesc = u"""Fitted to 37 distances.
""",
    longDesc = 
u"""
[<Entry index=153 label="OradH">, <Entry index=204 label="CsjCCC">]
[<Entry index=5 label="C_methane">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=1 label="H2">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=5 label="C_methane">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=204 label="CsjCCC">]
[<Entry index=153 label="OradH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=151 label="Cb_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=2 label="C_H">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=71 label="C_methyl">, <Entry index=204 label="CsjCCC">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=204 label="CsjCCC">]
[<Entry index=12 label="CsOHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=1 label="H2">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=7 label="CsCHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=204 label="CsjCCC">]
[<Entry index=2 label="C_H">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=161 label="OOH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=4 label="Csnorad_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=204 label="CsjCCC">]
[<Entry index=2 label="C_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=161 label="OOH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=1 label="H2">, <Entry index=204 label="CsjCCC">]
[<Entry index=12 label="CsOHHH">, <Entry index=226 label="Csj/Cs/Cs/O">]
""",
)

entry(
    index = 204,
    label = "CsjCCC",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    C  ux {1,S}
3    C  ux {1,S}
4    C  ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.023803, 'd13': 0.035775, 'd23': 0.00794},
        uncertainties = {'d12': 0.195386, 'd13': 0.181311, 'd23': 0.127588},
    ),
    shortDesc = u"""Fitted to 30 distances.
""",
    longDesc = 
u"""
[<Entry index=8 label="C/H3/Cs">, <Entry index=204 label="CsjCCC">]
[<Entry index=2 label="C_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=153 label="OradH">, <Entry index=204 label="CsjCCC">]
[<Entry index=153 label="OradH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=1 label="H2">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=161 label="OOH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=1 label="H2">, <Entry index=204 label="CsjCCC">]
[<Entry index=2 label="C_H">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=151 label="Cb_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=4 label="Csnorad_H">, <Entry index=204 label="CsjCCC">]
[<Entry index=71 label="C_methyl">, <Entry index=204 label="CsjCCC">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=204 label="CsjCCC">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=204 label="CsjCCC">]
[<Entry index=7 label="CsCHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=12 label="CsOHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=204 label="CsjCCC">]
[<Entry index=5 label="C_methane">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
""",
)

entry(
    index = 205,
    label = "Csj/Cs/Cs/Cs",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cs ux {1,S}
3    Cs ux {1,S}
4    Cs ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.019634, 'd13': 0.024895, 'd23': -0.002815},
        uncertainties = {'d12': 0.142934, 'd13': 0.162789, 'd23': 0.129969},
    ),
    shortDesc = u"""Fitted to 13 distances.
""",
    longDesc = 
u"""
[<Entry index=7 label="CsCHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=153 label="OradH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=161 label="OOH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=5 label="C_methane">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=2 label="C_H">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=1 label="H2">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
[<Entry index=12 label="CsOHHH">, <Entry index=205 label="Csj/Cs/Cs/Cs">]
""",
)

entry(
    index = 206,
    label = "Csj/Cs/Cs/Cd",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cs ux {1,S}
3    Cs ux {1,S}
4    Cd ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 207,
    label = "Csj/Cs/Cs/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cs ux {1,S}
3    Cs ux {1,S}
4    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 208,
    label = "Csj/Cs/Cs/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cs ux {1,S}
3    Cs ux {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 209,
    label = "Csj/Cs/Cd/Cd",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cs ux {1,S}
3    Cd ux {1,S}
4    Cd ux {1,S}
""",
    kinetics = DistanceData(distances={}),
    shortDesc = u"""Fitted to 1 distances.
""",
)

entry(
    index = 210,
    label = "Csj/Cs/Cd/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cs ux {1,S}
3    Cd ux {1,S}
4    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 211,
    label = "Csj/Cs/Cd/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cs ux {1,S}
3    Cd ux {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 212,
    label = "Csj/Cs/Ct/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cs ux {1,S}
3    Ct u0 {1,S}
4    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 213,
    label = "Csj/Cs/Ct/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cs ux {1,S}
3    Ct u0 {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 214,
    label = "Csj/Cs/Cb/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cs ux {1,S}
3    Cb u0 {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 215,
    label = "Csj/Cd/Cd/Cd",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cd ux {1,S}
3    Cd ux {1,S}
4    Cd ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 216,
    label = "Csj/Cd/Cd/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cd ux {1,S}
3    Cd ux {1,S}
4    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 217,
    label = "Csj/Cd/Cd/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cd ux {1,S}
3    Cd ux {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 218,
    label = "Csj/Cd/Ct/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cd ux {1,S}
3    Ct u0 {1,S}
4    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 219,
    label = "Csj/Cd/Ct/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cd ux {1,S}
3    Ct u0 {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 220,
    label = "Csj/Cd/Cb/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cd ux {1,S}
3    Cb u0 {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 221,
    label = "Csj/Ct/Ct/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Ct u0 {1,S}
3    Ct u0 {1,S}
4    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 222,
    label = "Csj/Ct/Ct/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Ct u0 {1,S}
3    Ct u0 {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 223,
    label = "Csj/Ct/Cb/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Ct u0 {1,S}
3    Cb u0 {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 224,
    label = "Csj/Cb/Cb/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cb u0 {1,S}
3    Cb u0 {1,S}
4    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 225,
    label = "CsjCCO",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    C  ux {1,S}
3    C  ux {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.041638, 'd13': 0.063305, 'd23': 0.021459},
        uncertainties = {'d12': 0.147058, 'd13': 0.163779, 'd23': 0.104978},
    ),
    shortDesc = u"""Fitted to 7 distances.
""",
    longDesc = 
u"""
[<Entry index=1 label="H2">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=2 label="C_H">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=161 label="OOH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=12 label="CsOHHH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=5 label="C_methane">, <Entry index=226 label="Csj/Cs/Cs/O">]
""",
)

entry(
    index = 226,
    label = "Csj/Cs/Cs/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cs ux {1,S}
3    Cs ux {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.041638, 'd13': 0.063305, 'd23': 0.021459},
        uncertainties = {'d12': 0.147058, 'd13': 0.163779, 'd23': 0.104978},
    ),
    shortDesc = u"""Fitted to 7 distances.
""",
    longDesc = 
u"""
[<Entry index=1 label="H2">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=2 label="C_H">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=161 label="OOH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=12 label="CsOHHH">, <Entry index=226 label="Csj/Cs/Cs/O">]
[<Entry index=5 label="C_methane">, <Entry index=226 label="Csj/Cs/Cs/O">]
""",
)

entry(
    index = 227,
    label = "Csj/Cs/Cd/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cs ux {1,S}
3    Cd ux {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 228,
    label = "Csj/Cs/Ct/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cs ux {1,S}
3    Ct u0 {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 229,
    label = "Csj/Cs/Cb/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cs ux {1,S}
3    Cb u0 {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 230,
    label = "Csj/Cd/Cd/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cd ux {1,S}
3    Cd ux {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 231,
    label = "Csj/Cd/Ct/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cd ux {1,S}
3    Ct u0 {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 232,
    label = "Csj/Cd/Cb/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cd ux {1,S}
3    Cb u0 {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 233,
    label = "Csj/Ct/Ct/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Ct u0 {1,S}
3    Ct u0 {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 234,
    label = "Csj/Ct/Cb/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Ct u0 {1,S}
3    Cb u0 {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 235,
    label = "Csj/Cb/Cb/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cb u0 {1,S}
3    Cb u0 {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 236,
    label = "CsjCOO",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    C  ux {1,S}
3    O  ux {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 237,
    label = "Csj/Cs/O/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cs ux {1,S}
3    O  ux {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 238,
    label = "Csj/Cd/O/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cd ux {1,S}
3    O  ux {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 239,
    label = "Csj/Ct/O/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Ct u0 {1,S}
3    O  ux {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 240,
    label = "Csj/Cb/O/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    Cb u0 {1,S}
3    O  ux {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 241,
    label = "CsjOOO",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    O  ux {1,S}
3    O  ux {1,S}
4    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 242,
    label = "Cdj",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u1
""",
    kinetics = DistanceData(
        distances = {'d12': -0.053236, 'd13': 0.058406, 'd23': 0.109644},
        uncertainties = {'d12': 0.130945, 'd13': 0.141613, 'd23': 0.164317},
    ),
    shortDesc = u"""Fitted to 299 distances.
""",
    longDesc = 
u"""
[<Entry index=9 label="C/H3/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=5 label="C_methane">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=153 label="OradH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=151 label="Cb_H">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=150 label="Ct_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=161 label="OOH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=5 label="C_methane">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=1 label="H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=155 label="OHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=7 label="CsCHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=12 label="CsOHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=5 label="C_methane">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=151 label="Cb_H">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=151 label="Cb_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=153 label="OradH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=2 label="C_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=4 label="Csnorad_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=155 label="OHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=5 label="C_methane">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=1 label="H2">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=7 label="CsCHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=1 label="H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=5 label="C_methane">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=151 label="Cb_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=4 label="Csnorad_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=155 label="OHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=7 label="CsCHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=155 label="OHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=7 label="CsCHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=2 label="C_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=153 label="OradH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=151 label="Cb_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=180 label="Csj/Cs/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=1 label="H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=7 label="CsCHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=161 label="OOH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=153 label="OradH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=1 label="H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=12 label="CsOHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=161 label="OOH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=250 label="Cdj_CdsCt">]
""",
)

entry(
    index = 243,
    label = "Cdj_CR",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u1 {2,D}
2    C  ux {1,D}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.053236, 'd13': 0.058406, 'd23': 0.109644},
        uncertainties = {'d12': 0.130945, 'd13': 0.141613, 'd23': 0.164317},
    ),
    shortDesc = u"""Fitted to 299 distances.
""",
    longDesc = 
u"""
[<Entry index=9 label="C/H3/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=5 label="C_methane">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=153 label="OradH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=151 label="Cb_H">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=150 label="Ct_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=161 label="OOH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=5 label="C_methane">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=1 label="H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=155 label="OHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=7 label="CsCHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=12 label="CsOHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=5 label="C_methane">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=151 label="Cb_H">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=151 label="Cb_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=153 label="OradH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=2 label="C_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=4 label="Csnorad_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=155 label="OHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=5 label="C_methane">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=1 label="H2">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=7 label="CsCHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=1 label="H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=5 label="C_methane">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=151 label="Cb_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=4 label="Csnorad_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=155 label="OHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=7 label="CsCHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=155 label="OHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=7 label="CsCHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=2 label="C_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=153 label="OradH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=151 label="Cb_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=180 label="Csj/Cs/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=1 label="H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=7 label="CsCHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=161 label="OOH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=153 label="OradH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=1 label="H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=12 label="CsOHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=161 label="OOH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=250 label="Cdj_CdsCt">]
""",
)

entry(
    index = 244,
    label = "Cdj_CH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u1 {2,D} {3,S}
2    C  ux {1,D}
3    H  u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.056394, 'd13': 0.055348, 'd23': 0.110922},
        uncertainties = {'d12': 0.136214, 'd13': 0.152488, 'd23': 0.172713},
    ),
    shortDesc = u"""Fitted to 194 distances.
""",
    longDesc = 
u"""
[<Entry index=161 label="OOH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=5 label="C_methane">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=1 label="H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=155 label="OHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=151 label="Cb_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=2 label="C_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=7 label="CsCHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=4 label="Csnorad_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=155 label="OHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=7 label="CsCHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=153 label="OradH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=150 label="Ct_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=151 label="Cb_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=180 label="Csj/Cs/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=5 label="C_methane">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=1 label="H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=161 label="OOH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=245 label="Cdj_CdsH">]
""",
)

entry(
    index = 245,
    label = "Cdj_CdsH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u1 {2,D} {3,S}
2    Cd ux {1,D}
3    H  u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.079715, 'd13': 0.059965, 'd23': 0.136924},
        uncertainties = {'d12': 0.142745, 'd13': 0.151594, 'd23': 0.165827},
    ),
    shortDesc = u"""Fitted to 121 distances.
""",
    longDesc = 
u"""
[<Entry index=5 label="C_methane">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=1 label="H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=161 label="OOH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=12 label="CsOHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=155 label="OHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=153 label="OradH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=150 label="Ct_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=151 label="Cb_H">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=7 label="CsCHHH">, <Entry index=245 label="Cdj_CdsH">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=245 label="Cdj_CdsH">]
""",
)

entry(
    index = 246,
    label = "Cdj_CddH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd  u1 {2,D} {3,S}
2    Cdd u0 {1,D}
3    H   u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.018433, 'd13': 0.047832, 'd23': 0.068596},
        uncertainties = {'d12': 0.127989, 'd13': 0.157628, 'd23': 0.187679},
    ),
    shortDesc = u"""Fitted to 73 distances.
""",
    longDesc = 
u"""
[<Entry index=7 label="CsCHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=2 label="C_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=155 label="OHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=5 label="C_methane">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=1 label="H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=148 label="Cdrad_Cdd/H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=4 label="Csnorad_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=151 label="Cb_H">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=12 label="CsOHHH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=161 label="OOH">, <Entry index=246 label="Cdj_CddH">]
[<Entry index=180 label="Csj/Cs/H2">, <Entry index=246 label="Cdj_CddH">]
""",
)

entry(
    index = 247,
    label = "Cdj_CC",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u1 {2,D} {3,S}
2    C  ux {1,D}
3    C  ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.047209, 'd13': 0.064244, 'd23': 0.107204},
        uncertainties = {'d12': 0.122802, 'd13': 0.121351, 'd23': 0.150321},
    ),
    shortDesc = u"""Fitted to 105 distances.
""",
    longDesc = 
u"""
[<Entry index=9 label="C/H3/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=151 label="Cb_H">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=1 label="H2">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=7 label="CsCHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=7 label="CsCHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=5 label="C_methane">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=151 label="Cb_H">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=155 label="OHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=153 label="OradH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=4 label="Csnorad_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=155 label="OHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=1 label="H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=5 label="C_methane">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=151 label="Cb_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=2 label="C_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=153 label="OradH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=12 label="CsOHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=7 label="CsCHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=161 label="OOH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=153 label="OradH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=1 label="H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=12 label="CsOHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=5 label="C_methane">, <Entry index=250 label="Cdj_CdsCt">]
""",
)

entry(
    index = 248,
    label = "Cdj_CdsCs",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u1 {2,D} {3,S}
2    Cd ux {1,D}
3    Cs ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.036709, 'd13': 0.066328, 'd23': 0.098647},
        uncertainties = {'d12': 0.099156, 'd13': 0.135816, 'd23': 0.155782},
    ),
    shortDesc = u"""Fitted to 39 distances.
""",
    longDesc = 
u"""
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=5 label="C_methane">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=151 label="Cb_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=12 label="CsOHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=2 label="C_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=4 label="Csnorad_H">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=7 label="CsCHHH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=161 label="OOH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=153 label="OradH">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=1 label="H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=248 label="Cdj_CdsCs">]
[<Entry index=155 label="OHH">, <Entry index=248 label="Cdj_CdsCs">]
""",
)

entry(
    index = 249,
    label = "Cdj_CdsCd",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u1 {2,D} {3,S}
2    Cd ux {1,D}
3    Cd ux {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.058528, 'd13': 0.059481, 'd23': 0.113781},
        uncertainties = {'d12': 0.149419, 'd13': 0.12718, 'd23': 0.163402},
    ),
    shortDesc = u"""Fitted to 51 distances.
""",
    longDesc = 
u"""
[<Entry index=9 label="C/H3/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=1 label="H2">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=151 label="Cb_H">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=155 label="OHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=153 label="OradH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=12 label="CsOHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=7 label="CsCHHH">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=5 label="C_methane">, <Entry index=249 label="Cdj_CdsCd">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=249 label="Cdj_CdsCd">]
""",
)

entry(
    index = 250,
    label = "Cdj_CdsCt",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u1 {2,D} {3,S}
2    Cd ux {1,D}
3    Ct u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.028868, 'd13': 0.077539, 'd23': 0.102413},
        uncertainties = {'d12': 0.100849, 'd13': 0.076804, 'd23': 0.117},
    ),
    shortDesc = u"""Fitted to 15 distances.
""",
    longDesc = 
u"""
[<Entry index=151 label="Cb_H">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=1 label="H2">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=7 label="CsCHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=153 label="OradH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=12 label="CsOHHH">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=120 label="Cdnorad_H">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=250 label="Cdj_CdsCt">]
[<Entry index=5 label="C_methane">, <Entry index=250 label="Cdj_CdsCt">]
""",
)

entry(
    index = 251,
    label = "Cdj_CdsCb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u1 {2,D} {3,S}
2    Cd ux {1,D}
3    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 252,
    label = "Cdj_CddCs",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd  u1 {2,D} {3,S}
2    Cdd u0 {1,D}
3    Cs  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 253,
    label = "Cdj_CddCd",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd  u1 {2,D} {3,S}
2    Cdd u0 {1,D}
3    Cd  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 254,
    label = "Cdj_CddCt",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd  u1 {2,D} {3,S}
2    Cdd u0 {1,D}
3    Ct  u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 255,
    label = "Cdj_CddCb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd  u1 {2,D} {3,S}
2    Cdd u0 {1,D}
3    Cb  u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 256,
    label = "Cdj_CO",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u1 {2,D} {3,S}
2    C  ux {1,D}
3    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 257,
    label = "Cdj_CdsO",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u1 {2,D} {3,S}
2    Cd ux {1,D}
3    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 258,
    label = "Cdj_CddO",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd  u1 {2,D} {3,S}
2    Cdd u0 {1,D}
3    O   ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 259,
    label = "Cdj_OR",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u1 {2,D}
2    O  u0 {1,D}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 260,
    label = "Cdj_OH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u1 {2,D} {3,S}
2    O  u0 {1,D}
3    H  u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 261,
    label = "Cdj_OC",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u1 {2,D} {3,S}
2    O  u0 {1,D}
3    C  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 262,
    label = "Cdj_OCs",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u1 {2,D} {3,S}
2    O  u0 {1,D}
3    Cs ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 263,
    label = "Cdj_OCd",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u1 {2,D} {3,S}
2    O  u0 {1,D}
3    Cd ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 264,
    label = "Cdj_OCt",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u1 {2,D} {3,S}
2    O  u0 {1,D}
3    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 265,
    label = "Cdj_OCb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u1 {2,D} {3,S}
2    O  u0 {1,D}
3    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 266,
    label = "Cdj_OO",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u1 {2,D} {3,S}
2    O  u0 {1,D}
3    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 267,
    label = "Ctj",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Ct u1
""",
    kinetics = DistanceData(
        distances = {'d12': -0.229192, 'd13': 0.225081, 'd23': 0.484929},
        uncertainties = {'d12': 0.361609, 'd13': 0.314334, 'd23': 0.210854},
    ),
    shortDesc = u"""Fitted to 20 distances.
""",
    longDesc = 
u"""
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=268 label="CtjC">]
[<Entry index=150 label="Ct_H">, <Entry index=268 label="CtjC">]
[<Entry index=7 label="CsCHHH">, <Entry index=268 label="CtjC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=268 label="CtjC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=268 label="CtjC">]
""",
)

entry(
    index = 268,
    label = "CtjC",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Ct u1 {2,T}
2    C  ux {1,T}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.229192, 'd13': 0.225081, 'd23': 0.484929},
        uncertainties = {'d12': 0.361609, 'd13': 0.314334, 'd23': 0.210854},
    ),
    shortDesc = u"""Fitted to 20 distances.
""",
    longDesc = 
u"""
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=268 label="CtjC">]
[<Entry index=150 label="Ct_H">, <Entry index=268 label="CtjC">]
[<Entry index=7 label="CsCHHH">, <Entry index=268 label="CtjC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=268 label="CtjC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=268 label="CtjC">]
""",
)

entry(
    index = 269,
    label = "Cbj",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cb u1
""",
    kinetics = DistanceData(
        distances = {'d12': -0.071008, 'd13': 0.063079, 'd23': 0.129934},
        uncertainties = {'d12': 0.235533, 'd13': 0.180956, 'd23': 0.200428},
    ),
    shortDesc = u"""Fitted to 30 distances.
""",
    longDesc = 
u"""
[<Entry index=12 label="CsOHHH">, <Entry index=269 label="Cbj">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=269 label="Cbj">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=269 label="Cbj">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=269 label="Cbj">]
[<Entry index=150 label="Ct_H">, <Entry index=269 label="Cbj">]
[<Entry index=151 label="Cb_H">, <Entry index=269 label="Cbj">]
[<Entry index=76 label="Csrad/H/Ct/H">, <Entry index=269 label="Cbj">]
[<Entry index=155 label="OHH">, <Entry index=269 label="Cbj">]
[<Entry index=2 label="C_H">, <Entry index=269 label="Cbj">]
[<Entry index=7 label="CsCHHH">, <Entry index=269 label="Cbj">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=269 label="Cbj">]
""",
)

entry(
    index = 270,
    label = "Cjj",
    group = "OR{Csjj, Cdjj}",
    kinetics = DistanceData(
        distances = {'d12': 0.003825, 'd13': 0.035387, 'd23': 0.035612},
        uncertainties = {'d12': 0.136846, 'd13': 0.138615, 'd23': 0.131268},
    ),
    shortDesc = u"""Fitted to 58 distances.
""",
    longDesc = 
u"""
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=1 label="H2">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=1 label="H2">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=161 label="OOH">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=12 label="CsOHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=7 label="CsCHHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=153 label="OradH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=5 label="C_methane">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=155 label="OHH">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=5 label="C_methane">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=151 label="Cb_H">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=155 label="OHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=1 label="H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=150 label="Ct_H">, <Entry index=305 label="Cs_trip/Ct/H">]
""",
)

entry(
    index = 271,
    label = "Csjj",
    group = "OR{Cs_sing, Cs_trip}",
    kinetics = DistanceData(
        distances = {'d12': -0.008504, 'd13': 0.055988, 'd23': 0.057546},
        uncertainties = {'d12': 0.151802, 'd13': 0.112008, 'd23': 0.147575},
    ),
    shortDesc = u"""Fitted to 32 distances.
""",
    longDesc = 
u"""
[<Entry index=155 label="OHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=155 label="OHH">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=1 label="H2">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=161 label="OOH">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=1 label="H2">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=5 label="C_methane">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=150 label="Ct_H">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=151 label="Cb_H">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=153 label="OradH">, <Entry index=305 label="Cs_trip/Ct/H">]
""",
)

entry(
    index = 272,
    label = "Cs_sing",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 273,
    label = "Cs_singH2",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    H  u0 {1,S}
3    H  u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 274,
    label = "Cs_singRH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs  u0 p1 {2,S} {3,S}
2    H   u0 {1,S}
3    R!H ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 275,
    label = "Cs_singCH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    H  u0 {1,S}
3    C  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 276,
    label = "Cs_sing/Cs/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    H  u0 {1,S}
3    Cs ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 277,
    label = "Cs_sing/Cd/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    H  u0 {1,S}
3    Cd ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 278,
    label = "Cs_sing/Ct/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    H  u0 {1,S}
3    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 279,
    label = "Cs_sing/Cb/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    H  u0 {1,S}
3    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 280,
    label = "Cs_singOH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    H  u0 {1,S}
3    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 281,
    label = "Cs_singRR",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs  u0 p1 {2,S} {3,S}
2    R!H ux {1,S}
3    R!H ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 282,
    label = "Cs_singCC",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    C  ux {1,S}
3    C  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 283,
    label = "Cs_sing/Cs/Cs",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    Cs ux {1,S}
3    Cs ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 284,
    label = "Cs_sing/Cs/Cd",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    Cs ux {1,S}
3    Cd ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 285,
    label = "Cs_sing/Cs/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    Cs ux {1,S}
3    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 286,
    label = "Cs_sing/Cs/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    Cs ux {1,S}
3    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 287,
    label = "Cs_sing/Cd/Cd",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    Cd ux {1,S}
3    Cd ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 288,
    label = "Cs_sing/Cd/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    Cd ux {1,S}
3    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 289,
    label = "Cs_sing/Cd/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    Cd ux {1,S}
3    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 290,
    label = "Cs_sing/Ct/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    Ct u0 {1,S}
3    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 291,
    label = "Cs_sing/Ct/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    Ct u0 {1,S}
3    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 292,
    label = "Cs_sing/Cb/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    Cb u0 {1,S}
3    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 293,
    label = "Cs_singCO",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    C  ux {1,S}
3    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 294,
    label = "Cs_sing/Cs/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    Cs ux {1,S}
3    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 295,
    label = "Cs_sing/Cd/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    Cd ux {1,S}
3    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 296,
    label = "Cs_sing/Ct/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    Ct u0 {1,S}
3    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 297,
    label = "Cs_sing/Cb/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    Cb u0 {1,S}
3    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 298,
    label = "Cs_singOO",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u0 p1 {2,S} {3,S}
2    O  ux {1,S}
3    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 299,
    label = "Cs_trip",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0
""",
    kinetics = DistanceData(
        distances = {'d12': -0.008504, 'd13': 0.055988, 'd23': 0.057546},
        uncertainties = {'d12': 0.151802, 'd13': 0.112008, 'd23': 0.147575},
    ),
    shortDesc = u"""Fitted to 32 distances.
""",
    longDesc = 
u"""
[<Entry index=155 label="OHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=155 label="OHH">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=1 label="H2">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=7 label="CsCHHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=161 label="OOH">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=1 label="H2">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=5 label="C_methane">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=150 label="Ct_H">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=151 label="Cb_H">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=153 label="OradH">, <Entry index=305 label="Cs_trip/Ct/H">]
""",
)

entry(
    index = 300,
    label = "Cs_tripH2",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    H  u0 {1,S}
3    H  u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': -0.049004, 'd13': 0.064985, 'd23': 0.105767},
        uncertainties = {'d12': 0.092622, 'd13': 0.127194, 'd23': 0.145563},
    ),
    shortDesc = u"""Fitted to 17 distances.
""",
    longDesc = 
u"""
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=1 label="H2">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=10 label="C/H3/Ct">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=155 label="OHH">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=161 label="OOH">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=300 label="Cs_tripH2">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=300 label="Cs_tripH2">]
""",
)

entry(
    index = 301,
    label = "Cs_tripRH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs  u2 p0 {2,S} {3,S}
2    R!H ux {1,S}
3    H   u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.047572, 'd13': 0.04353, 'd23': -0.009221},
        uncertainties = {'d12': 0.215392, 'd13': 0.10861, 'd23': 0.168946},
    ),
    shortDesc = u"""Fitted to 15 distances.
""",
    longDesc = 
u"""
[<Entry index=155 label="OHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=153 label="OradH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=150 label="Ct_H">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=1 label="H2">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=5 label="C_methane">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=151 label="Cb_H">, <Entry index=305 label="Cs_trip/Ct/H">]
""",
)

entry(
    index = 302,
    label = "Cs_tripCH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    C  ux {1,S}
3    H  u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.047572, 'd13': 0.04353, 'd23': -0.009221},
        uncertainties = {'d12': 0.215392, 'd13': 0.10861, 'd23': 0.168946},
    ),
    shortDesc = u"""Fitted to 15 distances.
""",
    longDesc = 
u"""
[<Entry index=155 label="OHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=153 label="OradH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=150 label="Ct_H">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=1 label="H2">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=5 label="C_methane">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=151 label="Cb_H">, <Entry index=305 label="Cs_trip/Ct/H">]
""",
)

entry(
    index = 303,
    label = "Cs_trip/Cs/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    Cs ux {1,S}
3    H  u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 304,
    label = "Cs_trip/Cd/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    Cd ux {1,S}
3    H  u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 305,
    label = "Cs_trip/Ct/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    Ct u0 {1,S}
3    H  u0 {1,S}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.047572, 'd13': 0.04353, 'd23': -0.009221},
        uncertainties = {'d12': 0.215392, 'd13': 0.10861, 'd23': 0.168946},
    ),
    shortDesc = u"""Fitted to 15 distances.
""",
    longDesc = 
u"""
[<Entry index=155 label="OHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=153 label="OradH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=150 label="Ct_H">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=1 label="H2">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=5 label="C_methane">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=7 label="CsCHHH">, <Entry index=305 label="Cs_trip/Ct/H">]
[<Entry index=151 label="Cb_H">, <Entry index=305 label="Cs_trip/Ct/H">]
""",
)

entry(
    index = 306,
    label = "Cs_trip/Cb/H",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    Cb u0 {1,S}
3    H  u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 307,
    label = "Cs_tripOH",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    O  ux {1,S}
3    H  u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 308,
    label = "Cs_tripRR",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs  u2 p0 {2,S} {3,S}
2    R!H ux {1,S}
3    R!H ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 309,
    label = "Cs_tripCC",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    C  ux {1,S}
3    C  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 310,
    label = "Cs_trip/Cs/Cs",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    Cs ux {1,S}
3    Cs ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 311,
    label = "Cs_trip/Cs/Cd",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    Cs ux {1,S}
3    Cd ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 312,
    label = "Cs_trip/Cs/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    Cs ux {1,S}
3    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 313,
    label = "Cs_trip/Cs/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    Cs ux {1,S}
3    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 314,
    label = "Cs_trip/Cd/Cd",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    Cd ux {1,S}
3    Cd ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 315,
    label = "Cs_trip/Cd/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    Cd ux {1,S}
3    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 316,
    label = "Cs_trip/Cd/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    Cd ux {1,S}
3    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 317,
    label = "Cs_trip/Ct/Ct",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    Ct u0 {1,S}
3    Ct u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 318,
    label = "Cs_trip/Ct/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    Ct u0 {1,S}
3    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 319,
    label = "Cs_trip/Cb/Cb",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    Cb u0 {1,S}
3    Cb u0 {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 320,
    label = "Cs_tripCO",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    C  ux {1,S}
3    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 321,
    label = "Cs_trip/Cs/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    Cs ux {1,S}
3    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 322,
    label = "Cs_trip/Cd/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    Cd ux {1,S}
3    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 323,
    label = "Cs_trip/Ct/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    Ct ux {1,S}
3    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 324,
    label = "Cs_trip/Cb/O",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    Cb ux {1,S}
3    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 325,
    label = "Cs_tripOO",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cs u2 p0 {2,S} {3,S}
2    O  ux {1,S}
3    O  ux {1,S}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 326,
    label = "Cdjj",
    group = "OR{Cd_singletR, Cd_tripletR}",
    kinetics = DistanceData(
        distances = {'d12': 0.017082, 'd13': 0.013237, 'd23': 0.012031},
        uncertainties = {'d12': 0.126598, 'd13': 0.174576, 'd23': 0.11841},
    ),
    shortDesc = u"""Fitted to 26 distances.
""",
    longDesc = 
u"""
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=1 label="H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=12 label="CsOHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=7 label="CsCHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=5 label="C_methane">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=331 label="Cd_tripletC">]
""",
)

entry(
    index = 327,
    label = "Cd_singletR",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u0 p1
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 328,
    label = "Cd_singletC",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u0 p1 {2,D}
2    C  ux {1,D}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 329,
    label = "Cd_singletO",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u0 p1 {2,D}
2    O  ux {1,D}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 330,
    label = "Cd_tripletR",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u2 p0
""",
    kinetics = DistanceData(
        distances = {'d12': 0.017082, 'd13': 0.013237, 'd23': 0.012031},
        uncertainties = {'d12': 0.126598, 'd13': 0.174576, 'd23': 0.11841},
    ),
    shortDesc = u"""Fitted to 26 distances.
""",
    longDesc = 
u"""
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=1 label="H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=12 label="CsOHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=7 label="CsCHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=5 label="C_methane">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=331 label="Cd_tripletC">]
""",
)

entry(
    index = 331,
    label = "Cd_tripletC",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u2 p0 {2,D}
2    C  ux {1,D}
""",
    kinetics = DistanceData(
        distances = {'d12': 0.017082, 'd13': 0.013237, 'd23': 0.012031},
        uncertainties = {'d12': 0.126598, 'd13': 0.174576, 'd23': 0.11841},
    ),
    shortDesc = u"""Fitted to 26 distances.
""",
    longDesc = 
u"""
[<Entry index=16 label="C/H2/Cs/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=1 label="H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=12 label="CsOHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=7 label="CsCHHH">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=9 label="C/H3/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=5 label="C_methane">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=123 label="Cd_Cds/H2">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=19 label="C/H2/Cd/Cd">, <Entry index=331 label="Cd_tripletC">]
[<Entry index=124 label="Cd_Cdd/H2">, <Entry index=331 label="Cd_tripletC">]
""",
)

entry(
    index = 332,
    label = "Cd_tripletO",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 Cd u2 p0 {2,D}
2    O  ux {1,D}
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 333,
    label = "Cjjj",
    group = "OR{C_doubletR, C_quartetR}",
    kinetics = DistanceData(
        distances = {'d12': -0.120266, 'd13': -0.052196, 'd23': 0.086769},
        uncertainties = {'d12': 0.269497, 'd13': 0.445836, 'd23': 0.32099},
    ),
    shortDesc = u"""Fitted to 3 distances.
""",
    longDesc = 
u"""
[<Entry index=1 label="H2">, <Entry index=335 label="C_quartetR">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=335 label="C_quartetR">]
[<Entry index=155 label="OHH">, <Entry index=335 label="C_quartetR">]
""",
)

entry(
    index = 334,
    label = "C_doubletR",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 C u1 p1
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 335,
    label = "C_quartetR",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 C u3
""",
    kinetics = DistanceData(
        distances = {'d12': -0.120266, 'd13': -0.052196, 'd23': 0.086769},
        uncertainties = {'d12': 0.269497, 'd13': 0.445836, 'd23': 0.32099},
    ),
    shortDesc = u"""Fitted to 3 distances.
""",
    longDesc = 
u"""
[<Entry index=1 label="H2">, <Entry index=335 label="C_quartetR">]
[<Entry index=8 label="C/H3/Cs">, <Entry index=335 label="C_quartetR">]
[<Entry index=155 label="OHH">, <Entry index=335 label="C_quartetR">]
""",
)

entry(
    index = 336,
    label = "Cjjjj",
    group = "OR{C_quintet, C_triplet}",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 337,
    label = "C_quintet",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 C u4 p0
""",
    kinetics = DistanceData(distances={}),
)

entry(
    index = 338,
    label = "C_triplet",
    group = 
"""
multiplicity [1,2,3,4,5]
1 *3 C u2 p1
""",
    kinetics = DistanceData(distances={}),
)

tree(
"""
L1: X_H_or_Xrad_H_Xbirad_H_Xtrirad_H
    L2: H2
    L2: C_H
        L3: Cs_H
            L4: Csnorad_H
                L5: C_methane
                L5: CsRHHH
                    L6: CsCHHH
                        L7: C/H3/Cs
                        L7: C/H3/Cd
                        L7: C/H3/Ct
                        L7: C/H3/Cb
                    L6: CsOHHH
                L5: CsRRHH
                    L6: CsCCHH
                        L7: C/H2/Cs/Cs
                        L7: C/H2/Cs/Cd
                        L7: C/H2/Cs/Ct
                        L7: C/H2/Cs/Cb
                        L7: C/H2/Cd/Cd
                        L7: C/H2/Cd/Ct
                        L7: C/H2/Cd/Cb
                        L7: C/H2/Ct/Ct
                        L7: C/H2/Ct/Cb
                        L7: C/H2/Cb/Cb
                    L6: CsCOHH
                        L7: C/H2/Cs/O
                        L7: C/H2/Cd/O
                        L7: C/H2/Ct/O
                        L7: C/H2/Cb/O
                    L6: CsOOHH
                L5: CsRRRH
                    L6: CsCCCH
                        L7: C/H/Cs/Cs/Cs
                        L7: C/H/Cs/Cs/Cd
                        L7: C/H/Cs/Cs/Ct
                        L7: C/H/Cs/Cs/Cb
                        L7: C/H/Cs/Cd/Cd
                        L7: C/H/Cs/Cd/Ct
                        L7: C/H/Cs/Cd/Cb
                        L7: C/H/Cs/Ct/Ct
                        L7: C/H/Cs/Ct/Cb
                        L7: C/H/Cs/Cb/Cb
                        L7: C/H/Cd/Cd/Cd
                        L7: C/H/Cd/Cd/Ct
                        L7: C/H/Cd/Cd/Cb
                        L7: C/H/Cd/Ct/Ct
                        L7: C/H/Cd/Ct/Cb
                        L7: C/H/Cd/Cb/Cb
                        L7: C/H/Ct/Ct/Ct
                        L7: C/H/Ct/Ct/Cb
                        L7: C/H/Ct/Cb/Cb
                        L7: C/H/Cb/Cb/Cb
                    L6: CsCCOH
                        L7: C/H/Cs/Cs/O
                        L7: C/H/Cs/Cd/O
                        L7: C/H/Cs/Ct/O
                        L7: C/H/Cs/Cb/O
                        L7: C/H/Cd/Cd/O
                        L7: C/H/Cd/Ct/O
                        L7: C/H/Cd/Cb/O
                        L7: C/H/Ct/Ct/O
                        L7: C/H/Ct/Cb/O
                        L7: C/H/Cb/Cb/O
                    L6: CsCOOH
                        L7: C/H/Cs/O/O
                        L7: C/H/Cd/O/O
                        L7: C/H/Ct/O/O
                        L7: C/H/Cb/O/O
                    L6: CsOOOH
            L4: Csrad_H
                L5: C_methyl
                L5: CsradRH2
                    L6: CsradCHH
                        L7: Csrad/H/Cs/H
                        L7: Csrad/H/Cd/H
                        L7: Csrad/H/Ct/H
                        L7: Csrad/H/Cb/H
                    L6: CsradOH2
                L5: CsradRRH
                    L6: CsradCCH
                        L7: Csrad/Cs/Cs/H
                        L7: Csrad/Cs/Cd/H
                        L7: Csrad/Cs/Ct/H
                        L7: Csrad/Cs/Cb/H
                        L7: Csrad/Cd/Cd/H
                        L7: Csrad/Cd/Ct/H
                        L7: Csrad/Cd/Cb/H
                        L7: Csrad/Ct/Ct/H
                        L7: Csrad/Ct/Cb/H
                        L7: Csrad/Cb/Cb/H
                    L6: CsradCOH
                        L7: Csrad/Cs/O/H
                        L7: Csrad/Cd/O/H
                        L7: Csrad/Ct/O/H
                        L7: Csrad/Cb/O/H
                    L6: CsradOOH
            L4: CsbiradH
                L5: Cs_singletH
                    L6: Cs_singletHH
                    L6: Cs_singletRH
                        L7: C_singletCH
                            L8: C_singlet/Cs/H
                            L8: C_singlet/Cd/H
                            L8: C_singlet/Ct/H
                            L8: C_singlet/Cb/H
                        L7: C_singletOH
                L5: Cs_tripletH
                    L6: Cs_tripletHH
                    L6: Cs_tripletRH
                        L7: Cs_tripletCH
                            L8: C_triplet/Cs/H
                            L8: C_triplet/Cd/H
                            L8: C_triplet/Ct/H
                            L8: C_triplet/Cb/H
                        L7: Cs_tripletOH
            L4: CstriradH
                L5: Cdoublet_H
                L5: Cquartet_H
        L3: Cd_H
            L4: Cdnorad_H
                L5: Cd_C/R/H
                    L6: Cd_C/H2
                        L7: Cd_Cds/H2
                        L7: Cd_Cdd/H2
                    L6: Cd_C/C/H
                        L7: Cd_Cds/Cs/H
                        L7: Cd_Cds/Cd/H
                        L7: Cd_Cds/Ct/H
                        L7: Cd_Cds/Cb/H
                        L7: Cd_Cdd/Cs/H
                        L7: Cd_Cdd/Cd/H
                        L7: Cd_Cdd/Ct/H
                        L7: Cd_Cdd/Cb/H
                    L6: Cd_C/O/H
                        L7: Cd_Cds/O/H
                        L7: Cd_Cdd/O/H
                L5: Cd_O/R/H
                    L6: Cd_O/H2
                    L6: Cd_O/C/H
                        L7: Cd_O/Cs/H
                        L7: Cd_O/Cd/H
                        L7: Cd_O/Ct/H
                        L7: Cd_O/Cb/H
                    L6: Cd_O/O/H
            L4: Cdrad_H
                L5: Cdrad_C/H
                    L6: Cdrad_Cds/H
                    L6: Cdrad_Cdd/H
                L5: Cdrad_O/H
        L3: Ct_H
        L3: Cb_H
    L2: O_H
        L3: OradH
        L3: ORH
            L4: OHH
            L4: OCH
                L5: O/Cs/H
                L5: O/Cd/H
                L5: O/Ct/H
                L5: O/Cb/H
            L4: OOH
L1: Y_rad_birad_trirad_quadrad
    L2: Hrad
    L2: Orad
        L3: OjR
            L4: OjH
            L4: OjC
                L5: OjCs
                L5: OjCd
                L5: OjCt
                L5: OjCb
            L4: OjO
        L3: O_atom_triplet
    L2: Crad
        L3: Cj
            L4: Csj
                L5: Cs_methyl
                L5: CsjRH2
                    L6: CsjCH2
                        L7: Csj/Cs/H2
                        L7: Csj/Cd/H2
                        L7: Csj/Ct/H2
                        L7: Csj/Cb/H2
                    L6: CsjOH2
                L5: CsjRRH
                    L6: CsjCCH
                        L7: Csj/Cs/Cs/H
                        L7: Csj/Cs/Cd/H
                        L7: Csj/Cs/Ct/H
                        L7: Csj/Cs/Cb/H
                        L7: Csj/Cd/Cd/H
                        L7: Csj/Cd/Ct/H
                        L7: Csj/Cd/Cb/H
                        L7: Csj/Ct/Ct/H
                        L7: Csj/Ct/Cb/H
                        L7: Csj/Cb/Cb/H
                    L6: CsjCOH
                        L7: Csj/Cs/O/H
                        L7: Csj/Cd/O/H
                        L7: Csj/Ct/O/H
                        L7: Csj/Cb/O/H
                    L6: CsjOOH
                L5: CsjRRR
                    L6: CsjCCC
                        L7: Csj/Cs/Cs/Cs
                        L7: Csj/Cs/Cs/Cd
                        L7: Csj/Cs/Cs/Ct
                        L7: Csj/Cs/Cs/Cb
                        L7: Csj/Cs/Cd/Cd
                        L7: Csj/Cs/Cd/Ct
                        L7: Csj/Cs/Cd/Cb
                        L7: Csj/Cs/Ct/Ct
                        L7: Csj/Cs/Ct/Cb
                        L7: Csj/Cs/Cb/Cb
                        L7: Csj/Cd/Cd/Cd
                        L7: Csj/Cd/Cd/Ct
                        L7: Csj/Cd/Cd/Cb
                        L7: Csj/Cd/Ct/Ct
                        L7: Csj/Cd/Ct/Cb
                        L7: Csj/Cd/Cb/Cb
                        L7: Csj/Ct/Ct/Ct
                        L7: Csj/Ct/Ct/Cb
                        L7: Csj/Ct/Cb/Cb
                        L7: Csj/Cb/Cb/Cb
                    L6: CsjCCO
                        L7: Csj/Cs/Cs/O
                        L7: Csj/Cs/Cd/O
                        L7: Csj/Cs/Ct/O
                        L7: Csj/Cs/Cb/O
                        L7: Csj/Cd/Cd/O
                        L7: Csj/Cd/Ct/O
                        L7: Csj/Cd/Cb/O
                        L7: Csj/Ct/Ct/O
                        L7: Csj/Ct/Cb/O
                        L7: Csj/Cb/Cb/O
                    L6: CsjCOO
                        L7: Csj/Cs/O/O
                        L7: Csj/Cd/O/O
                        L7: Csj/Ct/O/O
                        L7: Csj/Cb/O/O
                    L6: CsjOOO
            L4: Cdj
                L5: Cdj_CR
                    L6: Cdj_CH
                        L7: Cdj_CdsH
                        L7: Cdj_CddH
                    L6: Cdj_CC
                        L7: Cdj_CdsCs
                        L7: Cdj_CdsCd
                        L7: Cdj_CdsCt
                        L7: Cdj_CdsCb
                        L7: Cdj_CddCs
                        L7: Cdj_CddCd
                        L7: Cdj_CddCt
                        L7: Cdj_CddCb
                    L6: Cdj_CO
                        L7: Cdj_CdsO
                        L7: Cdj_CddO
                L5: Cdj_OR
                    L6: Cdj_OH
                    L6: Cdj_OC
                        L7: Cdj_OCs
                        L7: Cdj_OCd
                        L7: Cdj_OCt
                        L7: Cdj_OCb
                    L6: Cdj_OO
            L4: Ctj
                L5: CtjC
            L4: Cbj
        L3: Cjj
            L4: Csjj
                L5: Cs_sing
                    L6: Cs_singH2
                    L6: Cs_singRH
                        L7: Cs_singCH
                            L8: Cs_sing/Cs/H
                            L8: Cs_sing/Cd/H
                            L8: Cs_sing/Ct/H
                            L8: Cs_sing/Cb/H
                        L7: Cs_singOH
                    L6: Cs_singRR
                        L7: Cs_singCC
                            L8: Cs_sing/Cs/Cs
                            L8: Cs_sing/Cs/Cd
                            L8: Cs_sing/Cs/Ct
                            L8: Cs_sing/Cs/Cb
                            L8: Cs_sing/Cd/Cd
                            L8: Cs_sing/Cd/Ct
                            L8: Cs_sing/Cd/Cb
                            L8: Cs_sing/Ct/Ct
                            L8: Cs_sing/Ct/Cb
                            L8: Cs_sing/Cb/Cb
                        L7: Cs_singCO
                            L8: Cs_sing/Cs/O
                            L8: Cs_sing/Cd/O
                            L8: Cs_sing/Ct/O
                            L8: Cs_sing/Cb/O
                        L7: Cs_singOO
                L5: Cs_trip
                    L6: Cs_tripH2
                    L6: Cs_tripRH
                        L7: Cs_tripCH
                            L8: Cs_trip/Cs/H
                            L8: Cs_trip/Cd/H
                            L8: Cs_trip/Ct/H
                            L8: Cs_trip/Cb/H
                        L7: Cs_tripOH
                    L6: Cs_tripRR
                        L7: Cs_tripCC
                            L8: Cs_trip/Cs/Cs
                            L8: Cs_trip/Cs/Cd
                            L8: Cs_trip/Cs/Ct
                            L8: Cs_trip/Cs/Cb
                            L8: Cs_trip/Cd/Cd
                            L8: Cs_trip/Cd/Ct
                            L8: Cs_trip/Cd/Cb
                            L8: Cs_trip/Ct/Ct
                            L8: Cs_trip/Ct/Cb
                            L8: Cs_trip/Cb/Cb
                        L7: Cs_tripCO
                            L8: Cs_trip/Cs/O
                            L8: Cs_trip/Cd/O
                            L8: Cs_trip/Ct/O
                            L8: Cs_trip/Cb/O
                        L7: Cs_tripOO
            L4: Cdjj
                L5: Cd_singletR
                    L6: Cd_singletC
                    L6: Cd_singletO
                L5: Cd_tripletR
                    L6: Cd_tripletC
                    L6: Cd_tripletO
        L3: Cjjj
            L4: C_doubletR
            L4: C_quartetR
        L3: Cjjjj
            L4: C_quintet
            L4: C_triplet
"""
)

