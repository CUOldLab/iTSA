from matplotlib_venn import *
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def RemoveNan(this_set):
#Use the fact that nan != nan to remove them from the overlap set
    this_set = (x for x in this_set if x==x)
    return set(this_set)

def GenerateVenns(df,suffix):
#Get the overlap of optimum-temperature iTSA and our lab's TPP
    interior = len(RemoveNan(set(df["iTSA_52"]).intersection(df["Webb_TPP"])))

    left = len(RemoveNan(set(df["iTSA_52"]) - set(df["Webb_TPP"])))
    right = len(RemoveNan(set(df["Webb_TPP"]) - set(df["iTSA_52"])))

    out = venn2(subsets = (left,right,interior),
                set_labels = (r"iTSA 52$^{\circ}$C "+"\n"+str(len(RemoveNan(df["iTSA_52"]))),"Webb TPP\n "+str(len(RemoveNan(df["Webb_TPP"])))))
    c = venn2_circles((left,right,interior),linestyle = 'solid')
    
    for text in out.set_labels:
        text.set_fontsize(18)
    for text in out.subset_labels:
        text.set_fontsize(18)

    plt.savefig("iTSA52_Webb_Venn"+suffix+".png")
    plt.savefig("iTSA52_Webb_Venn"+suffix+".svg")
    plt.show()


#Get the overlap of all-temperature iTSA and our lab's TPP
    interior = len(RemoveNan(set(df["iTSA_ALL"]).intersection(df["Webb_TPP"])))
        
    left = len(RemoveNan(set(df["iTSA_ALL"]) - set(df["Webb_TPP"])))
    right = len(RemoveNan(set(df["Webb_TPP"]) - set(df["iTSA_52"])))

    out = venn2(subsets=(left,right,interior),set_labels=(r"iTSA all"+"\n"+str(len(RemoveNan(df["iTSA_ALL"]))),"Webb TPP\n "+str(len(RemoveNan(df["Webb_TPP"])))))
    c = venn2_circles((left,right,interior),linestyle='solid')
        
    for text in out.set_labels:
        text.set_fontsize(18)
    for text in out.subset_labels:
        text.set_fontsize(18)
    plt.savefig("iTSAall_Webb_Venn"+suffix+".png")
    plt.savefig("iTSAall_Webb_Venn"+suffix+".svg")
    plt.show()


    #Get the overlap of our lab's TPP and Savitski's

    interior = len(RemoveNan(set(df["Savitski_TPP"]).intersection(df["Webb_TPP"])))

    left = len(RemoveNan(set(df["Savitski_TPP"]) - set(df["Webb_TPP"])))
    right = len(RemoveNan(set(df["Webb_TPP"]) - set(df["Savitski_TPP"])))

    out = venn2(subsets=(left,right,interior),set_labels=(r"Savitski TPP"+"\n"+str(len(RemoveNan(df["Savitski_TPP"]))),"Webb TPP\n "+str(len(RemoveNan(df["Webb_TPP"])))))
    c = venn2_circles((left,right,interior),linestyle='solid')

    for text in out.set_labels:
        text.set_fontsize(18)
    for text in out.subset_labels:
        text.set_fontsize(18)
    plt.savefig("Savitski_Webb_Venn"+suffix+".png")
    plt.savefig("Savitski_Webb_Venn"+suffix+".svg")
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
    
    out = venn3(subsets=subsets,set_labels=(r"iTSA 48$^{\circ}$C "+"\n"+str(len(RemoveNan(df["iTSA_48"]))),r"iTSA 52$^{\circ}$C "+"\n"+str(len(RemoveNan(df["iTSA_52"]))),r"iTSA 56$^{\circ}$C "+"\n"+str(len(RemoveNan(df["iTSA_56"])))))
    c = venn3_circles(subsets=subsets,linestyle='solid')
    for text in out.set_labels:
        text.set_fontsize(18)
    for text in out.subset_labels:
        try:
            text.set_fontsize(18)
        except:
            pass
    plt.savefig("iTSA_temps_3way_Venn"+suffix+".png")
    plt.savefig("iTSA_temps_3way_Venn"+suffix+".svg")
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

    out = venn3(subsets=subsets,set_labels=("iTSA All\n"+str(len(RemoveNan(df["iTSA_ALL"]))) ,"Kinobeads\n"+str(len(RemoveNan(df["Kinobeads"]))),"Webb TPP \n"+str(len(RemoveNan(df["Webb_TPP"])))))
    c = venn3_circles(subsets=subsets,linestyle='solid')
    for text in out.set_labels:
        text.set_fontsize(18)
    for text in out.subset_labels:
        try:
            text.set_fontsize(18)
        except:
            pass
    plt.savefig("iTSAall_Webb_kinobead_Venn"+suffix+".png")
    plt.savefig("iTSAall_Webb_kinobead_Venn"+suffix+".svg")
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

    out = venn3(subsets=subsets,set_labels=("Savitski TPP \n"+str(len(RemoveNan(df["Savitski_TPP"]))),"Kinobeads\n"+str(len(RemoveNan(df["Kinobeads"]))),"Webb TPP \n"+str(len(RemoveNan(df["Webb_TPP"])))))
    c = venn3_circles(subsets=subsets,linestyle='solid')
    
    for text in out.set_labels:
        text.set_fontsize(18)
    for text in out.subset_labels:
        try:
            text.set_fontsize(18)
        except:
            pass
    plt.savefig("Savitski_Webb_kinobead_Venn"+suffix+".png")
    plt.savefig("Savitski_Webb_kinobead_Venn"+suffix+".svg")
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
    
    out = venn3(subsets=subsets,set_labels=(r"iTSA 48$^{\circ}$C "+"\n"+str(len(RemoveNan(df["iTSA_48"]))),"Kinobeads\n"+str(len(RemoveNan(df["Kinobeads"]))),r"iTSA 52$^{\circ}$C "+"\n"+str(len(RemoveNan(df["iTSA_52"])))))
    c = venn3_circles(subsets=subsets,linestyle='solid')
    
    for text in out.set_labels:
        text.set_fontsize(18)
    for text in out.subset_labels:
        try:
            text.set_fontsize(18)
        except:
            pass
    plt.savefig("iTSA_48_52_kinobead_Venn"+suffix+".png")
    plt.savefig("iTSA_48_52_kinobead_Venn"+suffix+".svg")
    plt.show()

def main():
    df = pd.read_csv("iTSA_overlaps.csv",sep=",")
    GenerateVenns(df,suffix="")

    df = pd.read_csv("iTSA_Overlaps_kinase.csv",sep=",")
    GenerateVenns(df,suffix="_kinase")

    df = pd.read_csv("iTSA_Overlaps_nonkinase.csv",sep=",")
    GenerateVenns(df,suffix="_nonkinase")

if __name__ == "__main__":
    main()
