import os 
#ce code permet de renommer les fichiers fastq en leur nom de souche.
name_dict={"ERR12524192_1.fastq":"YJS8725_1.fastq","ERR12524192_2.fastq":"YJS8725_2.fastq","ERR12523370_1.fastq":"YJS7833_1.fastq","ERR12523370_2.fastq":"YJS7833_2.fastq","ERR12523500_1.fastq":"YJS7992_1.fastq","ERR12523500_2.fastq":"YJS7992_2.fastq",
            "ERR12524172_1.fastq":"YJS8705_1.fastq","ERR12524172_2.fastq":"YJS8705_2.fastq","ERR12523722_1.fastq":"YJS8245_1.fastq","ERR12523722_2.fastq":"YJS8245_2.fastq","ERR4801040_1.fastq":"YJS7888_1.fastq","ERR4801040_2.fastq":"YJS7888_2.fastq",
            "ERR12523624_1.fastq":"YJS8144_1.fastq","ERR12523624_2.fastq":"YJS8144_2.fastq","ERR12523954_1.fastq":"YJS8486_1.fastq","ERR12523954_2.fastq":"YJS8486_2.fastq","ERR12523793_1.fastq":"YJS8318_1.fastq","ERR12523793_2.fastq":"YJS8318_2.fastq",
            "ERR4801046_1.fastq":"YJS8060_1.fastq","ERR4801046_2.fastq":"YJS8060_2.fastq","ERR12523821_1.fastq":"YJS8349_1.fastq","ERR12523821_2.fastq":"YJS8349_2.fastq","ERR12523921_1.fastq":"YJS8453_1.fastq","ERR12523921_2.fastq":"YJS8453_2.fastq",
            "ERR12523593_1.fastq":"YJS8112_1.fastq","ERR12523593_2.fastq":"YJS8112_2.fastq","ERR12523434_1.fastq":"YJS7918_1.fastq","ERR12523434_2.fastq":"YJS7918_2.fastq","ERR12523876_1.fastq":"YJS8407_1.fastq","ERR12523876_2.fastq":"YJS8407_2.fastq",
            "ERR12523592_1.fastq":"YJS8111_1.fastq","ERR12523592_2.fastq":"YJS8111_2.fastq","ERR12524014_1.fastq":"YJS8546_1.fastq","ERR12524014_2.fastq":"YJS8546_2.fastq","ERR12523466_1.fastq":"YJS7955_1.fastq","ERR12523466_2.fastq":"YJS7955_2.fastq",
            "ERR12524017_1.fastq":"YJS8549_1.fastq","ERR12524017_2.fastq":"YJS8549_2.fastq","ERR12523437_1.fastq":"YJS7922_1.fastq","ERR12523437_2.fastq":"YJS7922_2.fastq","ERR12523885_1.fastq":"YJS8416_1.fastq","ERR12523885_2.fastq":"YJS8416_2.fastq",
            "ERR12524098_1.fastq":"YJS8630_1.fastq","ERR12524098_2.fastq":"YJS8630_2.fastq","ERR12523787_1.fastq":"YJS8312_1.fastq","ERR12523787_2.fastq":"YJS8312_2.fastq","ERR12523459_1.fastq":"YJS7947_1.fastq","ERR12523459_2.fastq":"YJS7947_2.fastq",
            "ERR12524101_1.fastq":"YJS8633_1.fastq","ERR12524101_2.fastq":"YJS8633_2.fastq","ERR12523438_1.fastq":"YJS7923_1.fastq","ERR12523438_2.fastq":"YJS7923_2.fastq","ERR12523619_1.fastq":"YJS8139_1.fastq","ERR12523619_2.fastq":"YJS8139_2.fastq",
            "ERR12523545_1.fastq":"YJS8046_1.fastq","ERR12523545_2.fastq":"YJS8046_2.fastq","ERR12524185_1.fastq":"YJS8718_1.fastq","ERR12524185_2.fastq":"YJS8718_2.fastq","ERR12524102_1.fastq":"YJS8634_1.fastq","ERR12524102_2.fastq":"YJS8634_2.fastq"}

folder = "/data/projet3/reads_paired_end/reads/"

for old_name, new_name in name_dict.items():
    old_path = os.path.join(folder, old_name)
    new_path = os.path.join(folder, new_name)
    if os.path.exists(old_path):
        os.rename(old_path, new_path)
        print(f" {old_name} → {new_name}")
    else:
        print(f" {old_name} introuvable")
