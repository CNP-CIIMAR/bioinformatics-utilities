import subprocess
import re


def Cluster(File, Threshold, Statistics):
    print("++++++++++> work on", File)
    ref_Hash = initialMatrix(File)
    print("Do blast...")
    subprocess.call(["formatdb", "-i", File, "-o", "T"])
    Call = "blastall -b 10000000 -p blastp -e 0.1 -d " + File + " -i " + File + " -o " + File + ".tmp"
    subprocess.call(Call, shell=True)
    print("Did blast of", File)
    print("Start the parsing...")
    ref_Hash = parseScore(ref_Hash, File + ".tmp")
    print("Parsing done")
    ref_Hash = balanceHash(ref_Hash)
    print("Do interpretaation")
    ref_Group = interpretateHash(ref_Hash)
    print("Print Result")
    writeResult(ref_Group, File)
    print("Done for", File)
    subprocess.call(["rm", File + ".tmp"])


def initialMatrix(File):
    ref_Hash = {}
    SequenceAmount, ref_Names = getAmountSeq(File)
    Amount = len(ref_Names)
    for i in range(Amount):
        for j in range(Amount):
            ref_Hash[ref_Names[i]][ref_Names[j]] = 0
    return ref_Hash


def balanceHash(ref_Hash):
    for k1 in ref_Hash.keys():
        for k2 in ref_Hash[k1].keys():
            if k1 in ref_Hash and k2 in ref_Hash[k1] and ref_Hash[k1][k2] != ref_Hash[k2][k1]:
                ref_Hash[k1][k2] = 1
                ref_Hash[k2][k1] = 1
                print("not symmatric...")
            elif ((k1 not in ref_Hash or k2 in ref_Hash[k1]) or
                  (k1 in ref_Hash and k2 not in ref_Hash[k1])):
                ref_Hash[k1][k2] = 1
                ref_Hash[k2][k1] = 1
                print("not symmatric...")
    print("Balanced Matrix")
    return ref_Hash


def parseBlast(ref_Hash, FileName):
    ref_Parsed = returnm8(FileName)
    for line in ref_Parsed:
        elements = line.split("\t")
        qName, sName, qlength, slength, alignmentlength, ident, sim, Evalue, ScoreBit = elements
        if float(ScoreBit) >= Threshold:
            ref_Hash[qName][sName] = 1
    print("Did the parsing on", FileName)
    return ref_Hash


def parseScore(ref_Hash, FileName):
    with open(FileName, "r") as f3:
        Subject = ""
        Query = ""
        Score = ""
        for line in f3:
            if re.match(r"^Query=\s(\S+)", line):
                Query = re.match(r"^Query=\s(\S+)", line).group(1)
            elif re.match(r"^>(\S+)", line):
                Subject = re.match(r"^>(\S+)", line).group(1)
            elif re.match(r"^\sScore\s=\s*\S+\sbits\s+\((\d+)\),\sExp", line):
                                Score = re.match(r"^\sScore\s=\s*\S+\sbits\s+\((\d+)\),\sExp", line).group(1)
                Score = int(Score)
                if Score > 0:
                    ref_Hash[Query][Subject] = Score
    print("Did the parsing on", FileName)
    return ref_Hash


def interpretateHash(ref_Hash):
    print("Start the interpretation...")
    AllGroup = {}
    counter = 0
    for k1 in ref_Hash.keys():
        k1_group = ""
        countergroup = 0
        for k2 in ref_Hash[k1].keys():
            if ref_Hash[k1][k2] >= 40:
                k1_group += " " + k2
                countergroup += 1
        if countergroup > 0:
            k1_group = str(countergroup) + k1_group
            if k1_group not in AllGroup.values():
                AllGroup[counter] = k1_group
                counter += 1
    ref_Group = {}
    print("Take the intersection of groups...")
    for k1 in AllGroup.keys():
        for k2 in AllGroup.keys():
            intersection = set(AllGroup[k1].split()) & set(AllGroup[k2].split())
            if len(intersection) > 0:
                AllGroup[k1] = AllGroup[k1] + " " + AllGroup[k2]
    counter = 0
    for k1 in AllGroup.keys():
        if k1 in AllGroup.keys():
            AllGroup[k1] = AllGroup[k1].split()
            if AllGroup[k1] not in ref_Group.values():
                ref_Group[counter] = AllGroup[k1]
                counter += 1
    print("Interpreted the matrix")
    return ref_Group


def getAmountSeq(File):
    with open(File, "r") as f:
        ref_Names = []
        for line in f:
            if re.match(r"^>", line):
                ref_Names.append(line.strip()[1:])
        SequenceAmount = len(ref_Names)
    print("Read the file ", File, " found ", SequenceAmount, " sequences")
    return SequenceAmount, ref_Names


def writeResult(ref_Group, File):
    with open(File + ".result", "w") as f2:
        for k1 in ref_Group.keys():
            f2.write(">Group" + str(k1) + "\n")
            for k2 in range(int(ref_Group[k1][0])):
                f2.write(ref_Group[k1][k2 + 1] + "\n")
    print("result printed.")


if __name__ == "__main__":
    filename = "your_input_file.fasta"  # Put the filename of your input FASTA file here
    threshold_value = 50  # Put your desired threshold value here
    stats = False  # Set this to True if you want additional statistics
    Cluster(filename, threshold_value, stats)
