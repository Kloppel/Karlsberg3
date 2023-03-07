# -*- coding: utf-8 -*-
from audioop import rms

import re
import numpy as np
from matplotlib import pyplot as plt

from kbp2 import kbp_tools, analyse_md_pkas

def calc_rmsd(pka_list, exp_pka_list, resname_filter=None, segname_filter=None):
    # resname_filter = ['ASP', 'GLU', 'CYS', 'TYR', 'LYS', 'ARG', 'CTE', 'NTE', 'HIS']
    # resname_filter = ['ASP', 'GLU', 'CYS', 'TYR', 'LYS', 'ARG', 'CTE', 'HIS']
    # resname_filter = ['ASP', 'CYS', 'TYR', 'LYS', 'ARG', 'CTE', 'HIS']
    # resname_filter = ['ASP', 'GLU']
    # resname_filter = ['ASP']
    # resname_filter = ['GLU']
    # resname_filter = ['HIS', 'LYS']
    # resname_filter = ['NTE']
    # resname_filter = ['CTE']
    # resname_filter = ['CTE', 'NTE']
    # resname_filter = ['TYR']
    # resname_filter = ['LYS']
    # resname_filter = ['HIS']

    if type(pka_list) is not list:
        pka_list = [pka_list]
    if type(exp_pka_list) is not list:
        exp_pka_list = [exp_pka_list]
    assert(len(pka_list) == len(exp_pka_list))

    use_deviation = False
    # To get deviation, instead of RMSD for individual structures. Is helpful, if there is only one measured
    # pKA available,since it directly show the shift/deviation.
    #if len(exp_pka_list) == 1:
    #    use_deviation = True


    diffs = []
    for pkas, exp_pkas in zip(pka_list, exp_pka_list):
        diff = {}
        for residue_descr_kbp in pkas.keys():
            residue_descr = kbp_tools.get_real_residue(residue_descr_kbp)
            # diff[residue_descr] = None
            if resname_filter is not None:
                resname_kbp = re.split(r'[-_]', residue_descr_kbp)[0]
                resname = re.split(r'[-_]', residue_descr)[0]
                if (resname_kbp not in resname_filter) and (resname not in resname_filter):
                    continue
            if segname_filter is not None:
                segname = re.split(r'[-_]', residue_descr_kbp)[2]
                if (segname not in segname_filter):
                    continue
            pka = pkas[residue_descr_kbp]
            if residue_descr in exp_pkas:
                exp_pka = exp_pkas[residue_descr]
                diff[residue_descr] = pka - exp_pka
        diffs.append(diff)

    rmsds = []
    total_rmsd = 0.0
    total_nr_of_pkas = 0
    for diff in diffs:
        rmsd = 0.0
        nr_of_pkas = 0
        for delta_pka in diff.itervalues():
            if not use_deviation:
                rmsd += delta_pka**2
                total_rmsd += delta_pka**2
            else:
                rmsd += delta_pka
                total_rmsd += delta_pka

            nr_of_pkas += 1
            total_nr_of_pkas += 1
        if nr_of_pkas != 0:
            if not use_deviation:
                rmsd = np.sqrt(rmsd / float(nr_of_pkas))
            else:
                rmsd = rmsd / float(nr_of_pkas)
        else:
            rmsd = 0.0
        rmsds.append(rmsd)
    if total_nr_of_pkas != 0:
        if not use_deviation:
            total_rmsd = np.sqrt(total_rmsd / float(total_nr_of_pkas))
        else:
            total_rmsd = total_rmsd / float(total_nr_of_pkas)
    else:
        total_rmsd = 0.0

    return total_rmsd, rmsds, diffs

def interactive_plot(x, y_list, names, title='Tension', show=True):
    """
    @param x: []
    @param y: [[]]
    @param names: [description for outer list of y]
    @param title: Title of plot
    @param show: call.plt.show()
    @return: None
    """
    nr_of_lines = len(y_list)
    assert(len(y_list) == len(names))

    colors = plt.cm.Paired(np.linspace(0, 1, nr_of_lines))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.ion()

    line_names = {}
    lines = []

    for i, y, name in zip(range(nr_of_lines), y_list, names):
        line, = ax.plot(x, y, 'x-', color=colors[i], picker=5)
        line_names[line] = name
        lines.append(line)
    leg = ax.legend(names, loc='best', prop={'size': 7}, ncol=4)
    leg.draw_frame(False)
    plt.title(title)
    ymin, ymax = plt.ylim()
    plt.ylim([ymin - 0.1, 1.0])

    lined = dict()
    for legline, origline in zip(leg.get_lines(), lines):
        legline.set_picker(5)  # 5 pts tolerance
        lined[legline] = origline

    onpick = OnPick(title, line_names, lined, lines, ax, plt, names)

    fig.canvas.mpl_connect('pick_event', onpick)

    if show:
        plt.show(True)

class OnPick(object):
    def __init__(self, title, line_names, lined, lines, ax, plt, residue_list_filtered):
        self.line_names = line_names
        self.lined = lined
        self.ax = ax
        self.plt = plt
        self.residue_list_filtered = residue_list_filtered
        self.title = title
        self.lines = lines
    def __call__(self, event):
        from matplotlib.artist import setp

        for line in self.line_names.keys():
            setp(line, linewidth=1, marker='x')

        # lined[legline]
        if self.line_names.has_key(event.artist):
            thisLine = event.artist
            residue = self.line_names[thisLine]
        elif self.lined.has_key(event.artist):
            thisLine = self.lined[event.artist]
            residue = self.line_names[thisLine]

        setp(thisLine, linewidth=2, marker='o')
        # event.artist.set_label('---')

        # leg = ax.legend(residues, loc='best', prop={'size': 7}, ncol=4)
        leg = self.ax.legend(self.residue_list_filtered, loc='best', prop={'size': 7}, ncol=4)
        leg.draw_frame(False)
        # self.lined = dict()
        for legline, origline in zip(leg.get_lines(), self.lines):
            legline.set_picker(5)  # 5 pts tolerance
            self.lined[legline] = origline

        if self.title is not None:
            self.plt.title(self.title + ' - ' +  residue)
        else:
            self.plt.title(residue)


def print_project_rmsds(project_results, verbose=True, plot=True, output_str=[], resname_filter=None,
                        max_def_str_list=None):
    """
    output_str should be an empty list, the results are written in there.
    """

    delta_exp = []
    delta_calc = []
    delta_calc_old = []
    delta_residue_descr = []
    delta_names = []

    if verbose:
        print("Total RMSDs:")
    all_pkas_list = []
    all_exp_pkas_list = []
    all_old_pkas_list = []
    nr_of_residues = 0
    project_names = []
    for project_result in project_results:
        project_name = project_result.name
        project_names.append(project_name)
        pkas = project_result.pkas
        exp_pkas = project_result.exp_pkas
        old_pkas = project_result.old_pkas

        total_rmsd_list = []
        rmsds_list = []
        diffs_list = []
        total_rmsd_old_list = []
        rmsds_old_list = []
        diffs_old_list = []
        if project_name in ['2zta_cm2', '2zta_cm3', '1hng_cm', '2zta']:

            for segname in ['A', 'B']:
                total_rmsd, rmsds, diffs = calc_rmsd(pkas, exp_pkas, resname_filter=resname_filter
                                                     , segname_filter=segname)
                total_rmsd_list.append(total_rmsd)
                rmsds_list.append(rmsds)
                diffs_list.append(diffs)
                if (old_pkas is not None) and (type(old_pkas) is not list):
                    total_rmsd_old, rmsds_old, diffs_old = calc_rmsd(old_pkas, exp_pkas, resname_filter=resname_filter,
                                                                     segname_filter=segname)
                    total_rmsd_old_list.append(total_rmsd_old)
                    rmsds_old_list.append(rmsds_old)
                    diffs_old_list.append(diffs_old)
            # Average pkas
            pkas_avg = {}
            old_pkas_avg = {}
            for residue_descr_kbp, pka in pkas.iteritems():
                resname, resid, segname = re.split(r'[-_]', residue_descr_kbp)
                if segname == 'B':
                    segname = 'A'
                residue_descr_kbp_combined = "%s-%s_%s" % (resname, resid, segname)
                if residue_descr_kbp_combined in pkas_avg:
                    pkas_avg[residue_descr_kbp_combined] += pka * 0.5
                    if (old_pkas is not None) and (type(old_pkas) is not list):
                        old_pkas_avg[residue_descr_kbp_combined] += old_pkas[residue_descr_kbp] * 0.5
                else:
                    pkas_avg[residue_descr_kbp_combined] = pka * 0.5
                    if (old_pkas is not None) and (type(old_pkas) is not list):
                        old_pkas_avg[residue_descr_kbp_combined] = old_pkas[residue_descr_kbp] * 0.5
            pkas = pkas_avg
            if (old_pkas is not None) and (type(old_pkas) is not list):
                old_pkas = old_pkas_avg
        else:
            total_rmsd, rmsds, diffs = calc_rmsd(pkas, exp_pkas, resname_filter=resname_filter)
            if (old_pkas is not None) and (type(old_pkas) is not list):
                total_rmsd_old, rmsds_old, diffs_old = calc_rmsd(old_pkas, exp_pkas, resname_filter=resname_filter)

        all_pkas_list.append(pkas)
        all_exp_pkas_list.append(exp_pkas)
        if (old_pkas is not None) and (type(old_pkas) is not list):
            all_old_pkas_list.append(old_pkas)
        # nr_of_residues += len(exp_pkas)
        nr_of_residues_in_project = len(diffs[0])
        nr_of_residues += nr_of_residues_in_project

        if total_rmsd_list:
            # print total_rmsd_list
            total_rmsd_str = "%5.2f/%.2f" % tuple(total_rmsd_list)
            if (old_pkas is not None) and (type(old_pkas) is not list):
                total_rmsd_old_str = "%.2f/%.2f" % tuple(total_rmsd_old_list)
        else:
            total_rmsd_str = "%5.2f" % total_rmsd
            if (old_pkas is not None) and (type(old_pkas) is not list):
                total_rmsd_old_str = "%.2f" % total_rmsd_old
        if (old_pkas is not None) and (type(old_pkas) is not list):
            output_str.append((project_name, nr_of_residues_in_project, total_rmsd_str, total_rmsd_old_str))
            if verbose:
                print("  %8s (%i) : %s (KBP+: %s)" % (project_name, nr_of_residues_in_project,  total_rmsd_str, total_rmsd_old_str))
        else:
            output_str.append((project_name, nr_of_residues_in_project, total_rmsd_str, ''))
            if verbose:
                print("  %8s (%i): %5.2f" % (project_name, nr_of_residues_in_project, total_rmsd))
        for residue_descr in pkas:
            pka_calc = pkas[residue_descr]

            if (old_pkas is not None) and (type(old_pkas) is not list):
                if residue_descr not in old_pkas:
                    resname, resid , segname = re.split('[-_]', residue_descr)
                    # To allow mismatch between KB+ and KB+MD for the dPHS SNase variants.
                    if int(resid) in [45, 46, 48, 49, 52]:
                        continue
                    print "Residue not found in old pKa set:"
                    print project_name
                    print residue_descr
                    print old_pkas
                pka_calc_old = old_pkas[residue_descr]
            real_residue_descr = kbp_tools.get_real_residue(residue_descr)
            if real_residue_descr in exp_pkas:
                exp_pka = exp_pkas[real_residue_descr]
                # titratable_residues = project_result.detailed_project_results.descr.titratable_residues

                resname = re.split(r'[-_]', residue_descr)[0]
                # if resname in ['EPP', 'DPP']:
                # if resname in ['LYS']:
                #     continue
                # pka_mod = titratable_residues[resname][2]['pka'] + pka_mod_offset
                # pka_mod = abs(pka_mod)
                # print residue_descr
                # print pka_mod

                delta_calc_value = pka_calc
                if (old_pkas is not None) and (type(old_pkas) is not list):
                    delta_calc_old_value = pka_calc_old
                delta_exp_value = exp_pka
                # delta_calc_value = pka_calc - pka_mod
                # delta_calc_old_value = pka_calc_old - pka_mod
                # delta_exp_value = exp_pka - pka_mod

                delta_calc.append(delta_calc_value)
                if (old_pkas is not None) and (type(old_pkas) is not list):
                    delta_calc_old.append(delta_calc_old_value)
                delta_exp.append(delta_exp_value)
                delta_residue_descr.append(real_residue_descr)
                delta_names.append(project_name)


    total_rmsd, rmsds, diffs = calc_rmsd(all_pkas_list, all_exp_pkas_list, resname_filter=resname_filter)
    if (old_pkas is not None) and (type(old_pkas) is not list):
        total_rmsd_old, rmsds_old, diffs_old = calc_rmsd(all_old_pkas_list, all_exp_pkas_list,
                                                         resname_filter=resname_filter)
    if verbose:
        print("SASA factor: %.2f" % project_results[0].sasa_factor)

        if (old_pkas is not None) and (type(old_pkas) is not list):
            print("Total RMSD of %i residues: %.2f (KBP+: %.2f)" % (nr_of_residues, total_rmsd, total_rmsd_old))
        else:
            print("Total RMSD of %i residues: %.2f" % (nr_of_residues, total_rmsd))

    max_diff = 0
    max_diff_residue = None
    max_name = None
    all_diffs = []
    for name, diff in zip(project_names, diffs):
        for residue, delta_pka in diff.iteritems():
            # print delta_pka
            all_diffs.append(abs(delta_pka))
            if abs(max_diff) < abs(delta_pka):
                max_diff_residue = residue
                max_diff = delta_pka
                max_name = name
    if max_def_str_list is not None:
        max_def_str = "%s__%s__%.2f" % (max_name, max_diff_residue, max_diff)
        max_def_str_list.append(max_def_str)
    if verbose:
        print("Maximum deviation: %s %s : %.2f" % (max_name, max_diff_residue, max_diff))

    if (old_pkas is not None) and (type(old_pkas) is not list):
        max_diff = 0
        max_diff_residue = None
        max_name = None
        all_diffs_old = []
        for name, diff in zip(project_names, diffs_old):
            for residue, delta_pka in diff.iteritems():
                all_diffs_old.append(abs(delta_pka))
                if abs(max_diff) < abs(delta_pka):
                    max_diff_residue = residue
                    max_diff = delta_pka
                    max_name = name
        if max_def_str_list is not None:
            max_def_str_old = "%s__%s__%.2f" % (max_name, max_diff_residue, max_diff)
            max_def_str_list.append(max_def_str_old)
        if verbose:
            print("Maximum old deviation: %s %s : %.2f" % (max_name, max_diff_residue, max_diff))
    else:
        if max_def_str_list is not None:
            max_def_str_list.append(None)


    # bins = []
    # c = 0.0
    # while c < 10:
    #     bins.append(c)
    #     c += 0.5

    # plt.figure()
    # plt.hist(all_diffs, bins=bins)
    # plt.title('new')
    # plt.figure()
    # plt.hist(all_diffs_old, bins=bins)
    # plt.title('old')

    if plot:
        plt.figure()
        min_exp_delta = min(delta_exp)
        max_exp_delta = max(delta_exp)
        line = [min_exp_delta, max_exp_delta]
        plt.plot(delta_exp, delta_calc, 'ob')
        plt.plot(delta_exp, delta_calc_old, 'or')
        plt.plot(line, line, '-b')
        plt.legend(['new', 'old'], loc='best')
        plt.xlabel('experimental pka')
        plt.ylabel('calculated pka')

    # print "exp,calc,calc_old"
    # for exp, calc, calc_old, residue_descr in zip(delta_exp, delta_calc, delta_calc_old, delta_residue_descr):
    #     print "%.2f,%.2f,%.2f,%s" % (exp, calc, calc_old, residue_descr)
    # print line
    # print
    # print

    # exp_minus_mod = {}
    # calc_minus_mod = {}
    # old_calc_minus_mod = {}
    # calc_minus_exp = {}
    # old_calc_minus_exp = {}
    # titratable_residues = project_result.detailed_project_results.descr.titratable_residues
    # legend_printed = False
    # for exp, calc, calc_old, residue_descr in zip(delta_exp, delta_calc, delta_calc_old, delta_residue_descr):
    #     if not legend_printed:
    #         print "delta exp, delta calc, delta calc_old"
    #         legend_printed = True
    #     resname = re.split(r'[-_]', residue_descr)[0]
    #     resname_kb = kbp_tools.get_kbp_resname(resname)
    #     if titratable_residues[resname_kb][1]['pka'] < -50:
    #         # print titratable_residues[resname_kb]
    #         pka_mod = abs(titratable_residues[resname_kb][2]['pka'] + 100)
    #     else:
    #         pka_mod = abs(titratable_residues[resname_kb][1]['pka'])
    #     print "%.2f,%.2f,%.2f,%s" % (exp-pka_mod, calc-pka_mod, calc_old-pka_mod, residue_descr)
    #     if resname not in exp_minus_mod:
    #         exp_minus_mod[resname] = [0, 0]
    #         calc_minus_mod[resname] = [0, 0]
    #         old_calc_minus_mod[resname] = [0, 0]
    #         calc_minus_exp[resname] = [0, 0]
    #         old_calc_minus_exp[resname] = [0, 0]
    #     exp_minus_mod[resname][0] += 1
    #     exp_minus_mod[resname][1] += exp-pka_mod
    #
    #     calc_minus_mod[resname][0] += 1
    #     calc_minus_mod[resname][1] += calc-pka_mod
    #     old_calc_minus_mod[resname][0] += 1
    #     old_calc_minus_mod[resname][1] += calc_old-pka_mod
    #
    #     calc_minus_exp[resname][0] += 1
    #     calc_minus_exp[resname][1] += calc-exp
    #     old_calc_minus_exp[resname][0] += 1
    #     old_calc_minus_exp[resname][1] += calc_old-exp
    #
    # print "OLD"
    # for resname in old_calc_minus_mod.keys():
    #     count, shift = old_calc_minus_mod[resname]
    #     print "%s : %f" % (resname, shift / count)
    # print "calc - exp"
    # for resname in old_calc_minus_exp.keys():
    #     count, shift = old_calc_minus_exp[resname]
    #     print "%s : %f" % (resname, shift / count)
    # print
    #
    # print "NEW"
    # print "exp - mod"
    # for resname in exp_minus_mod.keys():
    #     count, shift = exp_minus_mod[resname]
    #     print "%s : %f" % (resname, shift / count)
    # print "calc - mod"
    # for resname in calc_minus_mod.keys():
    #     count, shift = calc_minus_mod[resname]
    #     print "%s : %f" % (resname, shift / count)
    # print "calc - exp"
    # for resname in calc_minus_exp.keys():
    #     count, shift = calc_minus_exp[resname]
    #     print "%s : %f" % (resname, shift / count)
    # print
    # print

    # for exp, calc, calc_old, residue_descr, name in zip(delta_exp, delta_calc, delta_calc_old, delta_residue_descr, delta_names):
    #     if abs(exp - calc_old) > 2.5:
    #         print "%.2f" % (exp - calc_old, )
    #         print "%.2f,%.2f,%.2f,%s,%s" % (exp, calc, calc_old, residue_descr, name)

    if (old_pkas is None) or (type(old_pkas) is list):
        total_rmsd_old = None

    return total_rmsd, rmsds, diffs, total_rmsd_old

def print_total_rmsds(pkas_list, exp_pkas_list, descriptions, old_pkas=None, resname_filter=None, plot=False):
    delta_exp = []
    delta_calc = []
    delta_calc_old = []

    # resname_filter = ['ASP', 'GLU']
    # resname_filter = ['LYS']

    print("Total RMSDs:")
    all_pkas_list = []
    all_exp_pkas_list = []
    all_old_pkas_list = []
    nr_of_residues = 0

    if type(pkas_list) is not list:
        pkas_list = [pkas_list]
    if type(exp_pkas_list) is not list:
        exp_pkas_list = [exp_pkas_list]
    if type(descriptions) is not list:
        descriptions = [descriptions]
    for pkas, exp_pkas, description in zip(pkas_list, exp_pkas_list, descriptions):
        total_rmsd, rmsds, diffs = calc_rmsd(pkas, exp_pkas, resname_filter=resname_filter)
        if (old_pkas is not None) and (type(old_pkas) is not list):
            total_rmsd_old, rmsds_old, diffs_old = calc_rmsd(old_pkas, exp_pkas, resname_filter=resname_filter)

        all_pkas_list.append(pkas)
        all_exp_pkas_list.append(exp_pkas)
        if (old_pkas is not None) and (type(old_pkas) is not list):
            all_old_pkas_list.append(old_pkas)
        # nr_of_residues += len(exp_pkas)
        nr_of_residues += len(diffs[0])

        if (old_pkas is not None) and (type(old_pkas) is not list):
            print("  %8s : %5.2f (prev: %.2f)" % (description, total_rmsd, total_rmsd_old))
        else:
            print("  %8s : %5.2f" % (description, total_rmsd))

        for residue_descr in pkas:
            pka_calc = pkas[residue_descr]
            if (old_pkas is not None) and (type(old_pkas) is not list):
                pka_calc_old = old_pkas[residue_descr]
            real_residue_descr = kbp_tools.get_real_residue(residue_descr)
            if real_residue_descr in exp_pkas:
                exp_pka = exp_pkas[real_residue_descr]

                delta_calc_value = pka_calc
                if (old_pkas is not None) and (type(old_pkas) is not list):
                    delta_calc_old_value = pka_calc_old
                delta_exp_value = exp_pka

                delta_calc.append(delta_calc_value)
                if (old_pkas is not None) and (type(old_pkas) is not list):
                    delta_calc_old.append(delta_calc_old_value)
                delta_exp.append(delta_exp_value)


    total_rmsd, rmsds, diffs = calc_rmsd(all_pkas_list, all_exp_pkas_list, resname_filter=resname_filter)
    if (old_pkas is not None) and (type(old_pkas) is not list):
        total_rmsd_old, rmsds_old, diffs_old = calc_rmsd(all_old_pkas_list, all_exp_pkas_list,
                                                         resname_filter=resname_filter)
    if (old_pkas is not None) and (type(old_pkas) is not list):
        print("Total RMSD of %i residues: %.2f (prev: %.2f)" % (nr_of_residues, total_rmsd, total_rmsd_old))
    else:
        print("Total RMSD of %i residues: %.2f" % (nr_of_residues, total_rmsd))

    all_diffs = []
    for diff in diffs:
        for delta_pka in diff.itervalues():
            # print delta_pka
            all_diffs.append(abs(delta_pka))
    all_diffs_old = []
    # if (old_pkas is not None) and (type(old_pkas) is not list):
    #     for diff in diffs_old:
    #         for delta_pka in diff.itervalues():
    #             all_diffs_old.append(abs(delta_pka))

    if plot:
        plt.figure()
        min_exp_delta = min(delta_exp)
        max_exp_delta = max(delta_exp)
        line = [min_exp_delta, max_exp_delta]
        plt.plot(delta_exp, delta_calc, 'ob')
        plt.plot(delta_exp, delta_calc_old, 'or')
        plt.plot(line, line, '-b')
        plt.legend(['new', 'old'], loc='best')
        plt.xlabel('experimental pka')
        plt.ylabel('calculated pka')

    return total_rmsd, rmsds, diffs

def select_project(project_results, selected_project_name):
    for project_result in project_results:
        if project_result.name == selected_project_name:
            return project_result
    return None

def print_project_pkas(project_result, all=False):
    project_name = project_result.name
    pkas = project_result.pkas
    exp_pkas = project_result.exp_pkas

    total_rmsd, rmsds, diffs = calc_rmsd(pkas, exp_pkas)
    diffs = diffs[0]

    if all:
        sorted_residue_list = get_sorted_residue_list(pkas.keys())
        sorted_residue_list = [kbp_tools.get_real_residue(x) for x in sorted_residue_list]
    else:
        sorted_residue_list = get_sorted_residue_list(diffs.keys())
    print("Results for project: %s" % project_name)
    for residue_descr in sorted_residue_list:
        kbp_residue_descr = kbp_tools.get_kbp_residue(residue_descr)
        pka = pkas[kbp_residue_descr]
        if residue_descr in exp_pkas:
            exp_pka = exp_pkas[residue_descr]
            # diff = diffs[residue_descr]
            diff = pka - exp_pka
            print "      %9s : %6.2f (%.2f) -> %6.2f" % (residue_descr, pka, exp_pka, diff)
        else:
            print "      %9s : %6.2f" % (residue_descr, pka)

def print_pkas(pkas, exp_pkas={}, description=None, all=False):
    if not exp_pkas:
        all = True

    if all:
        sorted_residue_list = get_sorted_residue_list(pkas.keys())
        sorted_residue_list = [kbp_tools.get_real_residue(x) for x in sorted_residue_list]
    else:
        total_rmsd, rmsds, diffs = calc_rmsd(pkas, exp_pkas)
        diffs = diffs[0]
        sorted_residue_list = get_sorted_residue_list(diffs.keys())
    if description is not None:
        print("Results for project: %s" % description)
    for residue_descr in sorted_residue_list:
        kbp_residue_descr = kbp_tools.get_kbp_residue(residue_descr)
        pka = pkas[kbp_residue_descr]
        if residue_descr in exp_pkas:
            exp_pka = exp_pkas[residue_descr]
            diff = pka - exp_pka
            print "      %9s : %6.2f (%5.2f) -> %6.2f" % (residue_descr, pka, exp_pka, diff)
        else:
            print "      %9s : %6.2f" % (residue_descr, pka)

def get_sorted_residue_list(residue_list):
    sorted_residue_list = list(residue_list)
    # Sort by resname.
    sorted_residue_list = sorted(residue_list)
    # Sort by segid.
    sorted_residue_list = sorted(sorted_residue_list, key=lambda residue: int(re.split('[_-]', residue)[1]))
    # Sorty by segname.
    sorted_residue_list = sorted(sorted_residue_list, key=lambda residue: re.split('[_-]', residue)[2])
    return sorted_residue_list

def plot_md_weight(project_result):
    detailed_project_results = project_result.detailed_project_results
    weights = detailed_project_results.weights
    ph_values = detailed_project_results.descr.ph_values
    task_names = detailed_project_results.get_task_names()

    nr_of_mds = np.shape(weights)[1]
    colors = plt.cm.Paired(np.linspace(0, 1, nr_of_mds))

    plt.figure()
    # print type(weights)
    for i, color in enumerate(colors):
        weight = weights[:, i]
        plt.plot(ph_values, weight, 'x-', color=color, linewidth=2)
    plt.legend(task_names)
    plt.ylim([-0.01, 1.01])

def plot_titration_curve(project_result, residue_descr):
    residue_nr = project_result.detailed_project_results.descr.get_residue_index(residue_descr)

    deprot_curve = project_result.detailed_project_results.combined_results.deprot_curves[residue_nr]
    ph_values = project_result.detailed_project_results.descr.ph_values
    plt.figure()
    plt.plot(ph_values, deprot_curve, 'rx-')
    plt.title("Deprotonation curve for " + residue_descr)
    plt.ylim([-0.01, 1.01])



def parse_propka_pkas(filename):
    pkas = {}
    f = open(filename)
    for line in f:
        entries = line.split()
        if len(entries) >= 6 and entries[5] == '%':
            resname, resid, segname, pka = entries[:4]
            resname = resname.replace('N+', 'NTE').replace('C-', 'CTE')
            resname = resname.replace('ASP', 'DPP').replace('GLU', 'EPP').replace('HIS', 'HSP')
            segname = segname.replace('_', 'A')
            residue_descr = "%s-%s_%s" % (resname, resid, segname)
            pka = pka.strip('*')
            pkas[residue_descr] = float(pka)
    f.close()

    return pkas


# ASP  19 A   4.22*   60 %    1.64  449   0.23    0   -0.85 THR  22 A   -0.01 ASP  19 A   -0.81 CA   CA B
# ASP  19 A                                            0.28 ASP  21 A   -0.48 ASP  21 A   -0.20 ARG  35 A
# ASP  19 A



# sasa_factors.append(sasa_factor)
# sasa_rmsds.append(total_rmsd)

# sasa_factors = []
# sasa_rmsds = []
# for sasa_factor in range(0, 130):
# sasa_factor = -4.8
# while sasa_factor < 2.0:
#     sasa_factor += 0.2
# for sasa_factor in [-0.4]:


# plt.figure()
# plt.plot(sasa_factors, sasa_rmsds)
# plt.show()
# for sasa_factor, sasa_rmsd in zip(sasa_factors, sasa_rmsds):
#     print("%.2f,%.2f" % (sasa_factor, sasa_rmsd))

