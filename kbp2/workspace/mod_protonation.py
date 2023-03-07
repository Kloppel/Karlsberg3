# -*- coding: utf-8 -*-

import cPickle as pickle
import kbp2

# protonation_pickle_filename     = '/scratch/scratch/tmeyer/md_pka/runs/general2/1ert_a/prot_tasks.pickle'
# protonation_pickle_filename     = '/scratch/scratch/tmeyer/md_pka/runs/general2/1ert_a/prot_tasks_min7_wo_26p.pickle'
# new_protonation_pickle_filename = '/scratch/scratch/tmeyer/md_pka/runs/general2/1ert_a/prot_tasks.pickle'


# protonation_pickle_filename     = '/scratch/scratch/tmeyer/md_pka/runs/general2/1xnb_35sw/prot_tasks.pickle'
# new_protonation_pickle_filename = '/scratch/scratch/tmeyer/md_pka/runs/general2/1xnb_35sw/prot_tasks.pickle'
# new_protonation_pickle_filename = '/scratch/scratch/tmeyer/md_pka/runs/general2/1xnb_35sw/prot_tasks_j172p.pickle'
protonation_pickle_filename     = '/scratch/scratch/tmeyer/md_pka/runs/general2/3bdc/ph7/prot_task.pickle'

new_protonation_pickle_filename = protonation_pickle_filename + '_previous'

f = open(protonation_pickle_filename)
(names, protonations) = pickle.load(f)
f.close

# f = open(new_protonation_pickle_filename, 'w')
# pickle.dump((names, protonations), f)
# f.close()
# # print names
# # print protonations
# new_protonations = {}
# for residue_tuple, state in protonations.iteritems():
#     if residue_tuple[0] not in ['CTE', 'NTE']:
#         new_protonations[residue_tuple] = state
#
# f = open(protonation_pickle_filename, 'w')
# pickle.dump((names, new_protonations), f)
# f.close()





# for name, protonation in zip(names, protonations):
#     if name != 'ph-10':
#         continue
#     print name
#     for residue_descr, state in protonation.iteritems():
#         print "%s : %s" % (residue_descr, str(state))
#
# protonation_pickle_filename     = '/scratch/scratch/tmeyer/md_pka/runs/general2/3bdc/prot_tasks_propka.pickle'
#
# f = open(protonation_pickle_filename)
# (names1, protonations1) = pickle.load(f)
# f.close
#
# for name, protonation, c in zip(names1, protonations1, range(len(names1))):
#     if name != 'ph-10_1':
#         continue
#
#     # print protonations[c].keys()
#     # print protonation.keys()
#
#     print name
#     for residue_descr, state in protonation.iteritems():
#         # print "%s : %s" % (residue_descr, str(state))
#
#         if state != protonations[c][residue_descr]:
#             print "%s : %s" % (residue_descr, str(state))
#             print protonations[c][residue_descr]




# print names
# names = ['ph5', 'ph7', 'ph-10', 'ph11_1']

# protonations[0][('NTE', 1, 'A')] = None
# protonations[0][('CTE', 33, 'A')] = None
# protonations[0][('GLU', 33, 'A')] = None
# protonations[0][('NTE', 1, 'B')] = None
# protonations[0][('CTE', 33, 'B')] = None
# protonations[0][('GLU', 33, 'B')] = None
#
# protonations[1][('NTE', 1, 'A')] = None
# protonations[1][('CTE', 33, 'A')] = None
# protonations[1][('GLU', 33, 'A')] = None
# protonations[1][('NTE', 1, 'B')] = None
# protonations[1][('CTE', 33, 'B')] = None
# protonations[1][('GLU', 33, 'B')] = None



# print names
#
# ph7_a_prot = kbp2.charmm.Charmm_manager.copy_titr_residue_dict( protonations[1] )
# print ph7_a_prot[('ASP', 26, 'A')]
# ph7_a_prot[('ASP', 26, 'A')] = 1

# names.append('ph7_1_26p')
# protonations.append(ph7_a_prot)


# ph5_172p = kbp2.charmm.Charmm_manager.copy_titr_residue_dict( protonations[1] )
# print ph5_172p[('GLU', 172, 'A')]
# ph5_172p[('GLU', 172, 'A')] = None
#
# names.append('ph5_172p')
# protonations.append(ph5_172p)
#
# index_to_remove = names.index('ph5')
# names.pop(index_to_remove)
# protonations.pop(index_to_remove)
#
# for name, protonation in zip(names, protonations):
#     print name
#     for residue_descr, state in protonation.iteritems():
#         if 172 in residue_descr:
#             print "%s : %s" % (residue_descr, str(state))



# f = open(new_protonation_pickle_filename, 'w')
# pickle.dump((names, protonations), f)
# f.close()
