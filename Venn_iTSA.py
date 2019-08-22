from matplotlib_venn import *
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd

df = pd.read_csv("iTSA_overlaps.csv",sep=",")

def RemoveNan(this_set):
#Use the fact that nan != nan to remove them from the overlap set
    this_set = (x for x in this_set if x==x)
    return set(this_set)


#Get the overlap of optimum-temperature iTSA and our lab's TPP
interior = len(RemoveNan(set(df["iTSA_52"]).intersection(df["Webb_TPP"])))

left = len(RemoveNan(set(df["iTSA_52"]) - set(df["Webb_TPP"])))
right = len(RemoveNan(set(df["Webb_TPP"]) - set(df["iTSA_52"])))

out = venn2(subsets=(left,right,interior),set_labels=(r"iTSA 52$^{\circ}$C "+"\n"+"95","Webb TPP\n 51"))
c = venn2_circles((left,right,interior),linestyle='solid')

for text in out.set_labels:
    text.set_fontsize(18)
for text in out.subset_labels:
    text.set_fontsize(18)

plt.savefig("iTSA52_Webb_Venn.png")
plt.show()


#Get the overlap of all-temperature iTSA and our lab's TPP
interior = len(RemoveNan(set(df["iTSA_ALL"]).intersection(df["Webb_TPP"])))

left = len(RemoveNan(set(df["iTSA_ALL"]) - set(df["Webb_TPP"])))
right = len(RemoveNan(set(df["Webb_TPP"]) - set(df["iTSA_52"])))

out = venn2(subsets=(left,right,interior),set_labels=(r"iTSA all"+"\n"+"130","Webb TPP\n 51"))
c = venn2_circles((left,right,interior),linestyle='solid')

for text in out.set_labels:
    text.set_fontsize(18)
for text in out.subset_labels:
    text.set_fontsize(18)
plt.savefig("iTSAall_Webb_Venn.png")
plt.show()


#Get the overlap of our lab's TPP and Savitski's

interior = len(RemoveNan(set(df["Savitski_TPP"]).intersection(df["Webb_TPP"])))

left = len(RemoveNan(set(df["Savitski_TPP"]) - set(df["Webb_TPP"])))
right = len(RemoveNan(set(df["Webb_TPP"]) - set(df["Savitski_TPP"])))

out = venn2(subsets=(left,right,interior),set_labels=(r"Savitski TPP"+"\n"+"60","Webb TPP\n 51"))
c = venn2_circles((left,right,interior),linestyle='solid')

for text in out.set_labels:
    text.set_fontsize(18)
for text in out.subset_labels:
    text.set_fontsize(18)
plt.savefig("Savitski_Webb_Venn.png")
plt.savefig("Savitski_Webb_Venn.svg")
plt.show()

###########################################
#########3-way Venns
#subsets=(Abc,aBc,ABc,abC,AbC,aBC,ABC)
###########################################
#The overlaps between the 3 iTSA runs

only_48 = len(RemoveNan(set(df["iTSA_48"]) - (set(df["iTSA_52"]) | set(df["iTSA_56"]))))
only_52 = len(RemoveNan(set(df["iTSA_52"]) - (set(df["iTSA_48"]) | set(df["iTSA_56"]))))
only_56 = len(RemoveNan(set(df["iTSA_56"]) - (set(df["iTSA_48"]) | set(df["iTSA_52"]))))
not_56 = len(RemoveNan(set(df["iTSA_48"]).intersection(set(df["iTSA_52"]) - set(df["iTSA_56"]))))
not_52 = len(RemoveNan(set(df["iTSA_48"]).intersection(set(df["iTSA_56"]) - set(df["iTSA_52"]))))
not_48 = len(RemoveNan(set(df["iTSA_56"]).intersection(set(df["iTSA_52"]) - set(df["iTSA_48"]))))
interior = len(RemoveNan(set(df["iTSA_48"]).intersection(set(df["iTSA_52"]),set(df["iTSA_56"]))))

subsets=(only_48,only_52,not_56,only_56,not_52,not_48,interior)

out = venn3(subsets=subsets,set_labels=(r"iTSA 48$^{\circ}$C "+"\n"+"72",r"iTSA 52$^{\circ}$C "+"\n"+"95",r"iTSA 56$^{\circ}$C "+"\n"+"66"))
c = venn3_circles(subsets=subsets,linestyle='solid')
for text in out.set_labels:
    text.set_fontsize(18)
for text in out.subset_labels:
    text.set_fontsize(18)
plt.savefig("iTSA_temps_3way_Venn.png")
plt.show()


#The comparison between iTSA_All, kinobeads, and Webb_TPP
only_iTSA = len(RemoveNan(set(df["iTSA_ALL"]) - (set(df["Kinobeads"]) | set(df["Webb_TPP"]))))
only_Kino = len(RemoveNan(set(df["Kinobeads"]) - (set(df["iTSA_ALL"]) | set(df["Webb_TPP"]))))
only_Webb = len(RemoveNan(set(df["Webb_TPP"]) - (set(df["iTSA_ALL"]) | set(df["Kinobeads"]))))
not_Webb = len(RemoveNan(set(df["iTSA_ALL"]).intersection(set(df["Kinobeads"])) - set(df["Webb_TPP"])))
not_Kino = len(RemoveNan(set(df["iTSA_ALL"]).intersection(set(df["Webb_TPP"])) - set(df["Kinobeads"])))
not_iTSA = len(RemoveNan(set(df["Webb_TPP"]).intersection(set(df["Kinobeads"])) - set(df["iTSA_ALL"])))
interior = len(RemoveNan(set(df["iTSA_ALL"]).intersection(set(df["Kinobeads"]),set(df["Webb_TPP"]))))

subsets = (only_iTSA,only_Kino,not_Webb,only_Webb,not_Kino,not_iTSA,interior)

out = venn3(subsets=subsets,set_labels=("iTSA All\n 130","Kinobeads\n 130","Webb TPP \n 51"))
c = venn3_circles(subsets=subsets,linestyle='solid')
for text in out.set_labels:
    text.set_fontsize(18)
for text in out.subset_labels:
    text.set_fontsize(18)
plt.savefig("iTSAall_Webb_kinobead_Venn.png")
plt.show()

#The comparison between Savitski TPP, kinobeads, and Webb_TPP
only_Savitski = len(RemoveNan(set(df["Savitski_TPP"]) - (set(df["Kinobeads"]) | set(df["Webb_TPP"]))))
only_Kino = len(RemoveNan(set(df["Kinobeads"]) - (set(df["Savitski_TPP"]) | set(df["Webb_TPP"]))))
only_Webb = len(RemoveNan(set(df["Webb_TPP"]) - (set(df["Savitski_TPP"]) | set(df["Kinobeads"]))))
not_Webb = len(RemoveNan(set(df["Savitski_TPP"]).intersection(set(df["Kinobeads"])) - set(df["Webb_TPP"])))
not_Kino = len(RemoveNan(set(df["Savitski_TPP"]).intersection(set(df["Webb_TPP"])) - set(df["Kinobeads"])))
not_Savitski = len(RemoveNan(set(df["Webb_TPP"]).intersection(set(df["Kinobeads"])) - set(df["Savitski_TPP"])))
interior = len(RemoveNan(set(df["Savitski_TPP"]).intersection(set(df["Kinobeads"]),set(df["Webb_TPP"]))))

subsets = (only_Savitski,only_Kino,not_Webb,only_Webb,not_Kino,not_Savitski,interior)

out = venn3(subsets=subsets,set_labels=("Savitski TPP \n 60","Kinobeads\n 130","Webb TPP \n 51"))
c = venn3_circles(subsets=subsets,linestyle='solid')
for text in out.set_labels:
    text.set_fontsize(18)
for text in out.subset_labels:
    text.set_fontsize(18)
plt.savefig("Savitski_Webb_kinobead_Venn.png")
plt.savefig("Savitski_Webb_kinobead_Venn.svg")
plt.show()



#The comparison between Savitski TPP, kinobeads, and Webb_TPP
only_Savitski = len(RemoveNan(set(df["iTSA_48"]) - (set(df["Kinobeads"]) | set(df["iTSA_52"]))))
only_Kino = len(RemoveNan(set(df["Kinobeads"]) - (set(df["iTSA_48"]) | set(df["iTSA_52"]))))
only_Webb = len(RemoveNan(set(df["iTSA_52"]) - (set(df["iTSA_56"]) | set(df["Kinobeads"]))))
not_Webb = len(RemoveNan(set(df["iTSA_48"]).intersection(set(df["Kinobeads"])) - set(df["iTSA_52"])))
not_Kino = len(RemoveNan(set(df["iTSA_48"]).intersection(set(df["iTSA_52"])) - set(df["Kinobeads"])))
not_Savitski = len(RemoveNan(set(df["iTSA_52"]).intersection(set(df["Kinobeads"])) - set(df["iTSA_48"])))
interior = len(RemoveNan(set(df["iTSA_48"]).intersection(set(df["Kinobeads"]),set(df["iTSA_52"]))))

subsets = (only_Savitski,only_Kino,not_Webb,only_Webb,not_Kino,not_Savitski,interior)

out = venn3(subsets=subsets,set_labels=(r"iTSA 48$^{\circ}$C "+"\n"+"72","Kinobeads\n 130",r"iTSA 52$^{\circ}$C "+"\n"+"95"))
c = venn3_circles(subsets=subsets,linestyle='solid')
for text in out.set_labels:
    text.set_fontsize(18)
for text in out.subset_labels:
    text.set_fontsize(18)
plt.savefig("iTSA_48_52_kinobead_Venn.png")
plt.savefig("iTSA_48_52_kinobead_Venn.svg")
plt.show()

