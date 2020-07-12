#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 21:42:04 2020

@author: lwoo0005
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_palette(sns.color_palette("tab20", 20))
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
import csv
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib import cm
from matplotlib.gridspec import GridSpec
from collections import OrderedDict

output_path='/Users/lwoo0005/Documents/Laura_stuff/H_pylori_project/'

#STEP 1: Find the transformation rate
# This is a massive simplification of Tim's model, but because the data are
# approximately linear, we can assume the transformation rate is linear as well

x_neutrals=np.array([a*23.1 for a in [0,0,0,2,2,2,4,4,4,6,6,6,7,7,7]])
x_neutral=np.array([float(a*7*3.33) for a in [0,2,4,6,7]])
y_neutral=np.array([a*0.01 for a in [0,0,0,0.958,1.271,0.617,2.666,2.567,2.266,3.301,4.087,2.176,3.284, 4.834,4.117]])
#y_neutral2=np.array([a*0.01 for a in [0,0,0,0.323,0.093,0.077,1.098,2.075,0.859,2.884,2.712,2.301,3.40,2.311,2.685]])
neutral_avg=[np.mean(y_neutral[0:3]), np.mean(y_neutral[3:6]), np.mean(y_neutral[6:9]), np.mean(y_neutral[9:12]), np.mean(y_neutral[12:15])]
rep_1=(y_neutral[0::3])
rep_2=(y_neutral[1::3])
rep_3=(y_neutral[2::3])

#This line to change to independent reps
neutral_avg=rep_1

"""
z = np.polyfit(x_neutral,neutral_avg,1, full=True)
#z = np.polyfit(x_neutrals,y_neutral,1, full=True)
"""



x_neutral = x_neutral[:,np.newaxis]
ga,SSE,_,_=np.linalg.lstsq(x_neutral, neutral_avg)

ga=float(ga)
CI_slope=float(SSE*1.96)
ga_up=float(ga+CI_slope)
ga_down=float(ga-CI_slope)
if ga_down <= 0:
    ga_down = 0
lower_bound=[ga_down*a*23.1 for a in [0,2,4,6,7]]
upper_bound=[ga_up*a*23.1 for a in [0,2,4,6,7]]
"""
f, ax =plt.subplots(nrows=1, ncols=1, figsize=(13,6))
#ax.plot(x_neutrals, y_neutral, label="Neutral allele average (observed)", marker='o', linestyle="None", color='green')
ax.plot(x_neutral, neutral_avg, label="Neutral allele (observed)", marker='o', linestyle="None", color='green')
ax.plot(x_neutral,ga*x_neutral, 'r-')
ax.fill_between(np.array([a*23.1 for a in [0,2,4,6,7]]), lower_bound, upper_bound, color='green', alpha=0.5)

plt.show()
"""
#plt.savefig(output_path+"/Poster_pics/Neutral_allele_trajectory.png", dpi=600)
#plt.close()


"""
f, ax =plt.subplots(nrows=1, ncols=1, figsize=(13,6))
ax.plot(x_neutral, neutral_avg, label="Neutral allele average (observed)", marker='o', linestyle="None", color='green')
ax.plot(x_neutral,ga*x_neutral, 'r-')
ax.fill_between(np.array([a*23.1 for a in [0,2,4,6,7]]), lower_bound, upper_bound, color='green', alpha=0.5)

plt.show()
plt.savefig(output_path+"/Poster_pics/Neutral_allele_trajectory.png", dpi=600)
plt.close()
"""


#Perform sel_coef calculation for all 4 time points and derive a CI for each
#gene, for each population

# Open the appropriate table and make a dataframe from the columns that
# indicate each abx-free timepoint of each of the H populations
# Split the main data into the 3 populations
# Define a function which takes t, g, z, prec, p_at_t
# Input the appropriate t value into the function

excel_file='/Users/lwoo0005/Documents/Laura_stuff/H_pylori_project/June_HP_mutations_compilation_p0.xlsx'
sheet = 'Newest'

pos= 0
gene= 31
HGT= 33
tpH1= 2
tpH2= 3
tpH3= 4
tpH4= 5

def df_collect_reps(excel_file, sheet, tp1, tp2, tp3, tp4, HGT, pos, gene):
    df = pd.read_excel(excel_file, sheet, index_col=None, na_values= ['']).iloc[:,[tp1, tp2, tp3, tp4, HGT, pos, gene]]
    df = df[df['Special?'].isin(['HGT', 'RdxA', 'FrxA'])]
    return df

H1= df_collect_reps(excel_file, sheet, tpH1, tpH2, tpH3, tpH4, HGT, pos, gene)
H2= df_collect_reps(excel_file, sheet, tpH1+5, tpH2+5, tpH3+5, tpH4+5, HGT, pos, gene)
H3= df_collect_reps(excel_file, sheet, tpH1+10, tpH2+10, tpH3+10, tpH4+10, HGT, pos, gene)

H_list=list([H1, H2, H3])

g_true = float(ga)

g_list = [g_true, 0.01, 0.001, 0.0001, 0.00001, 0.000001]

z = 0
#t_list = [float(t) for t in list(H1.columns[0:4])]
t_list = [float(t) for t in list(H1.columns[3:4])]

# The numbers defined in the linear space can be determined by a formula
# the paramters of which are given by the user-input prec value
s = np.linspace(-1, 1, 50001)

func1 = lambda s : (((z*s)+g)*(exp((s+g)*t))-(g*(1-z)))/(((z*s)+g)*(exp((s+g)*t))+(s*(1-z)))

t_fake=100
fake_t_list=[10, 100, 500, 1000]

p_eq_read_limited = float(0.001)
s_eq_read_limited = float(g_true)/(-1*p_eq_read_limited)

#set to 1 in 1000000
"""
p_eq_drift_limited = float(0.000001)
s_eq_drift_limited = float(g_true)/(-1*p_eq_drift_limited)
"""
s_eq_fix = g_true

xmin=s_eq_read_limited-0.1
xmax=s_eq_fix+0.1

width1=abs(s_eq_read_limited-xmin)
width2=abs(s_eq_fix-s_eq_read_limited)
width3=abs(xmax-s_eq_fix)

g_shade_list=[int(item) for item in list(np.linspace(60, 200, len(g_list)))]
t_shade_list=[int(item) for item in list(np.linspace(60, 200, len(fake_t_list)))]
sample_color_list=[cm.Oranges(150), cm.Purples(150), cm.Greens(150)]

plt.close()
fig=plt.figure(figsize=(10,12))
gs=GridSpec(2,2) # 2 rows, 2 columns
"""
ax=fig.add_subplot(gs[0,0])
ax2=fig.add_subplot(gs[0,1])
ax3=fig.add_subplot(gs[1,0:2])
"""
ax=fig.add_subplot(gs[1,0])
ax2=fig.add_subplot(gs[1,1])
ax3=fig.add_subplot(gs[0,0:2])


for g, shade in zip(g_list, g_shade_list):
    t=t_fake
    vals=func1(s)
    ax.plot(s, vals, label="%.0E" % (g), color=cm.Greys(shade))
    ax.set_xlim((xmin,xmax))
    ax.set_ylim((-0.0005,1))
    params = {'mathtext.default': 'regular' }
    plt.rcParams.update(params)
    ax.set_xlabel('Effective selection coefficient $(s_{e})$', fontsize=14)
    ax.set_ylabel("HGT allele frequency (P)", fontsize=14)
    ax.hlines(0,xmin=xmin,xmax=xmax, colors='grey', linestyles='dashed', lw=0.8)
    #rec1 = Rectangle([xmin,-0.002],height=1,width=width1)
    #rec2 = Rectangle([s_eq_read_limited,-0.002],height=1,width=width2)
    #rec3 = Rectangle([s_eq_fix,-0.002],height=1,width=width3)
    #rec = PatchCollection([rec1,rec2,rec3],facecolor=['red','orange','yellow'],alpha=0.04,edgecolor='none')
    #ax.add_collection(rec)
bounds=[]
for g in [ga_up, ga_down]:
    t=t_fake
    if g==0:
        bound_vals=[0 for c in s]
    else:
        bound_vals=func1(s)
    bounds.append(bound_vals)
ax.fill_between(s, bounds[0], bounds[1], color=cm.YlOrRd(170), alpha=0.2)
ax.get_lines()[0].set_color(cm.YlOrRd(170))
ax.set_yscale('symlog', linthreshy=0.000001)
#ax.set_yscale('linear')
leg=ax.legend(prop={'size': 10}, fancybox=True, ncol=2, title='Invasion rate at generation 100', fontsize=10, columnspacing=0.6, handlelength=1, labelspacing=0.3, handletextpad=0.6)
leg.get_frame().set_alpha(0.6)
leg._legend_box.align = "left"
leg.get_texts()[0].set_text('%.2E $\pm$ 95%% CI' % (g_true))
#leg.get_texts()[-1].set_text('No invasion;\nmutation rate: %.0E' % (g_list[-1]))
ax.set_title('b)', loc='left', fontsize=14)

#########
"""
plt.close
fig=plt.figure(figsize=(15,20))
gs=GridSpec(2,2) # 2 rows, 2 columns
ax=fig.add_subplot(gs[0,0])
ax2=fig.add_subplot(gs[0,1])
ax3=fig.add_subplot(gs[1,0:2])
"""
#######

for t, shade in zip(fake_t_list, t_shade_list):
    g=g_true
    vals=func1(s)
    ax2.plot(s, vals, label="%i" % (t), color=cm.Greys(shade))
    ax2.set_xlim((xmin,xmax))
    ax2.set_ylim((-0.01,0.02))
    params = {'mathtext.default': 'regular' }
    plt.rcParams.update(params)
    ax2.set_xlabel('Effective selection coefficient $(s_{e})$', fontsize=14)
    #ax2.set_ylabel("HGT-derived allele\nfrequency (P)", fontsize=14)
    ax2.hlines(0,xmin=xmin,xmax=xmax, colors='grey', linestyles='dashed', lw=0.8)
    rec1 = Rectangle([xmin,ax2.get_ylim()[0]],height=(ax2.get_ylim()[1]-ax2.get_ylim()[0]),width=width1)
    rec2 = Rectangle([s_eq_read_limited,ax2.get_ylim()[0]],height=(ax2.get_ylim()[1]-ax2.get_ylim()[0]),width=width2)
    rec3 = Rectangle([s_eq_fix,ax2.get_ylim()[0]],height=(ax2.get_ylim()[1]-ax2.get_ylim()[0]),width=width3)
    rec = PatchCollection([rec1,rec2,rec3],facecolor=[cm.YlOrRd(150),cm.YlOrRd(80),cm.YlOrRd(40)],alpha=0.1,edgecolor='none')
    ax2.add_collection(rec)
#ax2.get_lines()[0].set_color("cyan")
#ax2.set_yscale('symlog', linthreshy=0.0001)
ax2.set_yscale('linear')
custom_lines=[Line2D([0],[0], marker='s', markersize=12, mfc=cm.YlOrRd(150), mec='k', mew=0.8, alpha=0.2, lw=0),
                Line2D([0],[0], marker='s', markersize=12, mfc=cm.YlOrRd(80), mec='k', mew=0.8, alpha=0.2, lw=0),
                Line2D([0],[0], marker='s', markersize=12, mfc=cm.YlOrRd(40), mec='k', mew=0.8, alpha=0.2, lw=0)]
sec_leg = ax2.legend(custom_lines, ['Undetectable at equilibrium', 'Undetectable < equilibrium frequency < 1', 'Fixation at equilibrium'], loc=3, prop={'size': 10}, fancybox=True, labelspacing=0.3, borderpad=0.4, fontsize=10, handletextpad=0.6)
leg=ax2.legend(prop={'size': 10}, fancybox=True, ncol=4, title='Generation at experimental transformation rate', loc=2, columnspacing=0.8, fontsize=10, handlelength=1, borderpad=0.4, handletextpad=0.6)
sec_leg.get_frame().set_alpha(0.6)
leg.get_frame().set_alpha(0.6)
leg._legend_box.align = "left"
ax2.add_artist(sec_leg)
ax2.set_title('c)', loc='left', fontsize=14)  
#plt.savefig(output_path+"/Poster_pics/S_P_relationship_t100_all_gs.png", dpi=600)

def closest(lst, K):  
     idx = (np.abs(lst - K)).argmin() 
     return lst[idx] 
 
"""    
    
H_num = 1
for H, sample_color in zip(H_list, sample_color_list):
    all_ts=[]
    for t in t_list:
        z= 0
        g= g_true
        vals=func1(s)
        sel_coefs=[]
        if H_num > 1:
            col_H = H_num-1
            col_of_i=str(t)+'.'+str(col_H)
        else:
            col_of_i=float(t)
        H_items=[float(p) for index, p in H.loc[:,[col_of_i]].iterrows()]
        H_pos=[int(position) for index, position in H.iloc[:,[5]].iterrows()]
        #H_genes=[(str(ge).split('\\')[0]).split('\t')[1] for index, ge in H.iloc[:,[6]].iterrows()]
        #H_genes=[str(ge) for index, ge in H.iloc[:,[6]].iterrows()]
        H_genes=[(ge['gene']).encode('utf-8') for index, ge in H.iloc[:,[6]].iterrows()]
        for p, position, gel in zip(H_items, H_pos, H_genes):
            close_p = closest(vals, p)
            for sel, val in zip(s, vals):
                if val == close_p:
                    s_final = sel
                    ax3.scatter(gel, s_final, marker='o', alpha=0.3, edgecolors=sample_color, color =sample_color, label='Population H%i' % H_num)
                    sel_coefs.append(s_final)
        all_ts.append(sel_coefs)
    with open(str(output_path)+"/sel_coefs_by_Tim_rep_%i.tsv" % H_num, 'w') as op:
        op_writer=csv.writer(op, delimiter='\t')
        #op_writer.writerow(['Position', 'Gene', 'Time_1', 'Time_2', 'Time_3', 'Time_4'])
        op_writer.writerow(['Position', 'Gene', 'Time_4'])
        #for pos, ges, t1, t2, t3, t4 in zip(H_pos, H_genes, all_ts[0], all_ts[1], all_ts[2], all_ts[3]):
            #op_writer.writerow([pos, ges, t1, t2, t3, t4])
        for pos, ges, t4 in zip(H_pos, H_genes, all_ts[0]):
            op_writer.writerow([pos, ges, t4])
    H_num+=1
params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)
ax3.set_ylabel('Effective selection coefficient $(s_{e})$', fontsize=13)
ax3.set_xlabel('Gene', fontsize=13)
ax3.set_ylim(s_eq_read_limited,0.2)
handles, labels = ax3.get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
ax3.legend(by_label.values(), by_label.keys(), prop={'size': 8})
ax3.set_title('c)', loc='left', fontsize=15)
plt.show()

"""


low_bounds=[]
high_bounds=[]
H_num = 1
for H, sample_color in zip(H_list, sample_color_list):
    all_ts=[]
    for t in t_list:
        z= 0    
        g= g_true
        vals=func1(s)
        sel_coefs=[]
        if H_num > 1:
            col_H = H_num-1
            col_of_i=str(t)+'.'+str(col_H)
        else:
            col_of_i=float(t)
        H_items=[float(p) for index, p in H.loc[:,[col_of_i]].iterrows()]
        H_pos=[int(position) for index, position in H.iloc[:,[5]].iterrows()]
        #H_genes=[(str(ge).split('\\')[0]).split('\t')[1] for index, ge in H.iloc[:,[6]].iterrows()]
        #H_genes=[str(ge) for index, ge in H.iloc[:,[6]].iterrows()]
        #H_genes=[unicode((ge['gene']).encode('utf-8')) for index, ge in H.iloc[:,[6]].iterrows()]
        #H_genes=H_genes=['hopG', 'ccdA', 'dppC', 'HPP12_0753', 'cheA', 'IpxA', 'HPP12_1495', 'mod-5 (1)', 'mod-5 (2)']
        H_genes=['hopG', 'ccdA', 'dppC', 'frxA', 'HPP12_0753', 'rdxA (1)', 'rdxA (2)', 'rdxA (3)', 'rdxA (4)', 'rdxA (5)', 'rdxA (6)', 'rdxA (7)', 'rdxA (8)', 'rdxA (9)', 'rdxA (10)', 'rdxA (11)', 'rdxA (12)', 'rdxA (13)', 'rdxA (14)','rdxA (15)', 'rdxA (16)','rdxA (17)','rdxA (18)','rdxA (19)','rdxA (20)','rdxA (21)','rdxA (22)','rdxA (23)','cheA', 'IpxA', 'HPP12_1495', 'mod-5 (1)', 'mod-5 (2)']
        for p, position, i in zip(H_items, H_pos, range(1,len(H_genes)+1)):
            close_p = closest(vals, p)
            for sel, val in zip(s, vals):
                if val == close_p:
                    s_final = sel
                    #ax3.scatter(i, s_final, marker='o', alpha=0.3, edgecolors=sample_color, color =sample_color, label='Population H%i' % H_num)
                    sel_coefs.append(s_final)
            g=ga_down
            low_vals=func1(s)
            close_p_low=closest(low_vals,p)
            if ga_down==0:
                close_p_low=0
            for sel, low_val_trans in zip(s, low_vals):
                if low_val_trans == close_p_low:
                    s_low_trans = sel
                    high_bounds.append(sel)
            g=ga_up
            high_vals=func1(s)
            close_p_high=closest(high_vals,p)
            for sel, high_val_trans in zip(s, high_vals):
                if high_val_trans == close_p_high:
                    s_high_trans = sel
                    low_bounds.append(sel)
            CI_width=0.4
            CI_box = Rectangle(xy=[i-(0.5*CI_width),s_high_trans], width=CI_width,height=s_low_trans-s_high_trans)
            CI_boxy = PatchCollection([CI_box],facecolor=sample_color,alpha=0.3,edgecolor='none')
            ax3.add_collection(CI_boxy)
            g=g_true           
        all_ts.append(sel_coefs)
    with open(str(output_path)+"/sel_coefs_by_Tim_rep_%i.tsv" % H_num, 'w') as op:
        op_writer=csv.writer(op, delimiter='\t')
        #op_writer.writerow(['Position', 'Gene', 'Time_1', 'Time_2', 'Time_3', 'Time_4'])
        op_writer.writerow(['Position', 'Gene', 'Time_4'])
        #for pos, ges, t1, t2, t3, t4 in zip(H_pos, H_genes, all_ts[0], all_ts[1], all_ts[2], all_ts[3]):
            #op_writer.writerow([pos, ges, t1, t2, t3, t4])
        for pos, ges, t4 in zip(H_pos, H_genes, all_ts[0]):
            op_writer.writerow([pos, ges, t4])
    H_num+=1
params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)
ax3.set_ylabel('Effective selection\ncoefficient $(s_{e})$', fontsize=14)
ax3.set_xlabel('Gene', fontsize=14)
ax3.set_ylim(min(low_bounds)-0.01,max(high_bounds)+0.01)
#ax3.set_ylim(s_eq_read_limited,0.2)
ax3.set_xlim(range(1,len(H_genes)+1)[0]-(2*CI_width), range(1,len(H_genes)+1)[-1]+(2*CI_width))
#ax3.set_xlim(0.5,9.5)
"""
handles, labels = ax3.get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
ax3.legend(by_label.values(), by_label.keys(), prop={'size': 8})
"""
custom_lines_ax3=[Line2D([0],[0], marker='s', markersize=12, mfc=cm.Oranges(150), mec='k', mew=0.8, alpha=0.2, lw=0),
                Line2D([0],[0], marker='s', markersize=12, mfc=cm.Purples(150), mec='k', mew=0.8, alpha=0.2, lw=0),
                Line2D([0],[0], marker='s', markersize=12, mfc=cm.Greens(150), mec='k', mew=0.8, alpha=0.2, lw=0)]
leg_ax3 = ax3.legend(custom_lines_ax3, ['Population HGT 1 95% CI', 'Population HGT 2 95% CI', 'Population HGT 3 95% CI'], loc='best', prop={'size': 10}, fancybox=True, labelspacing=0.3, borderpad=0.4, handletextpad=0.6)
ax3.add_artist(leg_ax3)
ax3.set_title('a)', loc='left', fontsize=14)
ax3.set_xticks(range(1,len(H_genes)+1))
ax3.set_xticklabels(H_genes, fontstyle='italic', rotation=75, fontsize=9)
ax3.hlines(0,xmin=ax3.get_xlim()[0],xmax=ax3.get_xlim()[1], colors='grey', linestyles='dashed', lw=0.8)
"""
plt.subplots_adjust(
top=0.96,
bottom=0.17,
left=0.12,
right=0.98,
hspace=0.345,
wspace=0.21
)
"""
plt.subplots_adjust(
top=0.96,
bottom=0.07,
left=0.12,
right=0.98,
hspace=0.55,
wspace=0.21
)

plt.show()

plt.savefig(output_path+"/Poster_pics/three_panel_tims_sel_coefs_fig_with_metR_alleles_scaled_down_rep1_newlayout.png", dpi=600)
plt.savefig(output_path+"/Poster_pics/three_panel_tims_sel_coefs_fig_with_metR_alleles_scaled_down_rep1_newlayout.pdf", dpi=600)
print "All done--thanks for being the best!"