#!/home/nanje/miniconda3/bin/python

def optionR():
    """It reads in pdb files"""
    import os.path
    
    global file
    file =input("Please enter the name of the PDB file: ")
    
    global path
    path= '../../python-mini-project-nanjalaruth/Data/'
    global validfile
    validfile=path+file
    
    if os.path.exists(validfile) == True:
        print(file,"is a valid file")
        
        global files_read
        files_read=[]
        if file in files_read:
            print("file already loaded")
        else:
            files_read.append(file)
               
    else:
        print("Invalid file selected: Please try again!!")
        optionR()
        
        
"""main dictionary; the amino acid codon is represented by the key while the value represents its 1 letter translation"""
maindict={
'ALA':"A",'CYS':"C",'ASP':"D",'GLU':"E",'PHE':"F",'GLY':"G",'HIS':"H",'ILE':"I",'LYS':"K",'LEU':"L",'MET':"M",'ASN':"N",'PRO':"P",'GLN':"Q",'ARG':"R",'SER':"S",'THR':"T",'VAL':"V",'TRP':"W",'TYR':"Y"}


def optionS():
    """Search option"""
    
    print("These are the files that have been read",files_read)
    import re
    raw_sequence=[]
    seff=[]
    option=input("please select a file to work with:")
    valid_file=path + option
    if option in valid_file:
        with open(valid_file)as opt:
            for line in opt:
                if line.startswith("SEQRES"):
                    list=(str.split(line[19:70]))
                    for res in list:
                        seff.append(res)
                        raw_sequence.append(maindict[res])
        sequence=str(''.join(raw_sequence).replace(' ',''))
        pattern=re.compile(r'[GYL][A][PFW][TLVGMT]')
        rhyme=pattern.finditer(sequence)
        lyst = []
        for m in rhyme:
            m.group()
            x = m.span()
            y= x[0]
            z = x[1]
            for i in range(y, z):
                lyst.append(i)
        sequenc= []
        for i in sequence:
            sequenc.append(i)
        for c in lyst:
            indx = str(sequenc[c])
            lwer=indx.lower()
            sequenc[c]=lwer
        print("".join(sequenc))
    else: 
        print("invalid choice")
        optionS()

        
        
"""packages that come in handy during alignment and writing out sequences in FASTA format"""        
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC



def optionC():
    """writes out coordinate section of the files"""
    of= input("please select an output file:")
    pth="../../python-mini-project-nanjalaruth/Results/"
    vfile=pth + of
    with open(valid_file,'r')as ab:
        with open(vfile,'w')as uc:
            for line in ab:
                if line.startswith("ATOM"):
                    op=(line[32:55] + '\n')
                    uc.writelines(op)
                if line.startswith("HETATM"):
                    p=(line[28:56] + '\n')
                    uc.writelines(p)
    def soptions():
        kk=input("please select '5' to write out all atom, '6' to write out backbone atom or '7' to write out Alpha carbon atoms or '8' to exit:")
        if kk !='8':
            if kk in ('5','6','7'):
                if kk =='5':
                    """It writes out coordinate file for all atoms"""
                    grc=input("please provide an output file:")
                    pth="../../python-mini-project-nanjalaruth/Results/"
                    vfileq=pth + grc
                    with open(valid_file,'r')as aa:
                        with open( vfileq,'w')as cc:
                            for line in aa:
                                if line.startswith("ATOM"):
                                    op=(line[32:55] + '\n')
                                    cc.writelines(op)
                    soptions()
                elif kk =='6':
                    """It writes out coordinate file for backbone atoms"""
                    grc1=input("please provide an output file")
                    pth=('../../python-mini-project-nanjalaruth/Results/')
                    vfileq1=pth+grc1
                    with open(valid_file,'r')as ab:
                        with open(vfileq1,'w')as p:
                            for line in ab:
                                if line.startswith("ATOM"):
                                    if "CA" in line:
                                        op=(line[32:55] + '\n')
                                        p.writelines(op)
                                    elif "C" in line:
                                        op=(line[32:55] + '\n')
                                        p.writelines(op)
                                    elif "N" in line:
                                        op=(line[32:55] + '\n')
                                        p.writelines(op)
                                    elif "O" in line:
                                        op=(line[32:55] + '\n')
                                        p.writelines(op)
                    soptions()
                elif kk =='7':
                    """It writes out coordinate file-alpha carbon atom"""
                    grc2=input("please provide an output file")
                    pth=('../../python-mini-project-nanjalaruth/Results/')
                    vfileq2=pth+grc2
                    with open(valid_file,'r')as ab:
                        with open(vfileq2,'w')as k:
                            for line in ab:
                                if line.startswith("ATOM"):
                                    if "CA" in line:
                                        op=(line[32:55] + '\n')
                                        k.writelines(op)
                    soptions()
                elif kk =='8':
                    sys.exit()
                else:
                    print("Invalid option! try again")
                    soptions()
    soptions()

    
    

def optionT():
    """It write out sequence in Fasta format"""
    print("These are the files that have been read",files_read)
    oon=input("please select a file to work with")
    valifile=path + oon
    if oon in valifile:
        def soptions2():
            pk=input("please select '1' to write out SEQRES sequence, '2' to write out coordinate sequence, '3' to write out Alignment sequence or '4' to quit:")
            if pk !='4':
                if pk in ('1','2','3'):
                    if pk =='1':
                        """It writes out SEQRES sequences in Fasta format"""
                        sef2=[]
                        rsequence2=[]
                        pthe=('../../python-mini-project-nanjalaruth/Results/')
                        kem2=input("please select an output file")
                        opfile=pthe+kem2
                        with open(valifile,'r')as op2:
                            with open(opfile, 'w')as ss:
                                for sline in op2:
                                    if sline.startswith("SEQRES"):
                                        line=(str.split(sline[19:]))
                                        for rs2 in line:
                                            sef2.append(rs2)
                                            rsequence2.append(maindict[rs2])
                                seqc2=str(''.join(rsequence2).replace(' ',''))
                                Sseq = SeqRecord(Seq(seqc2,IUPAC.protein),id="RN",description=f'PDB file in {file}')
                                ss.write(Sseq.format('fasta'))
                                print(f'file successfully written to: {kem2}')
                        soptions2()
                    elif pk =='2':
                        """It writes out coordinate sequences in Fasta format"""
                        y=''
                        pthe=('../../python-mini-project-nanjalaruth/Results/')
                        kem1=input("please select an output file")
                        opfile1=pthe+kem1
                        with open('../../python-mini-project-nanjalaruth/Data/1zni.pdb','r')as op1:
                            with open(opfile1,'w')as q: 
                                for kline in op1:
                                    if kline.startswith("ATOM"):
                                        line=(kline[16:20])
                                        y=y+line
                                cseq = SeqRecord(Seq(str(y),IUPAC.protein),id="RN",description=f'PDB file in {file}')
                                q.write(cseq.format('fasta'))
                                print(f'file successfully written to: {kem1}')       
                        soptions2()
                    elif pk =='3':
                        """It writes out Alignment sequence in Fasta format"""
                        sef=[]
                        rsequence=[]
                        pthe=('../../python-mini-project-nanjalaruth/Results/')
                        kem=input("please select an output file")
                        opfile=pthe+kem
                        with open(valifile,'r')as op:
                            for sline in op:
                                if sline.startswith("SEQRES"):
                                    line=(str.split(sline[19:70]))
                                    for rs in line:
                                        sef.append(rs)
                                        rsequence.append(maindict[rs])
                            seqc=str(''.join(rsequence).replace(' ',''))
                            seqc_list =[]
                            for i in seqc:
                                seqc_list.append(i)
                        with open(valifile,'r')as op:
                            with open(opfile, 'w')as cc:
                                rem = []
                                for line in op: 
                                    if "REMARK 465" in line:
                                        k = line[10:15]
                                        if k == "     ":
                                            p = line[23:]
                                            if p[:4] != "    ":
                                                num= int(p)
                                                rem.append(num)
                                for q in rem:
                                    seqc_list[q-1]="-X-"
                                p=("".join(seqc_list))
                                Sseq = SeqRecord(Seq(p,IUPAC.protein),id="RN",description=f'PDB file in {file}')
                                cc.write(Sseq.format('fasta'))
                                print(f'file successfully written to: {kem}')
                        soptions2()
                    elif pk =='4':
                        sys.exit()
                else:
                    print("Invalid option! try again")
                    soptions2()
        soptions2()
    else:
        print("Invalid option!")
        optionT()

    
        
def optionW():
    """Write out option"""
    print("These are the files that have been read",files_read)
    opt=input("please select a file to work with:")
    global valid_file
    valid_file=path + opt
    def miniopt():
        ket=input("please select 'C' to write out coordinate file,'T' to write out sequence files or 'Q' to Quit:")
        if ket.capitalize() !='Q':
            if ket.capitalize()in('C','T'):
                if ket.capitalize()=='C':
                    optionC()
                    miniopt()
                elif ket.capitalize()=='T':
                    optionT()
                    miniopt()
                elif ket.capitalize()=='Q':
                    sys.exit()
            else:
                print("Invalid option!")
                miniopt()
    miniopt()




def optionI():
    """Information option"""
    print("These are the files that have been read",files_read)
    option=input("please select a file to work with")
    global valfile
    valfile=path + option
    if option in valfile:
        def minioption():
            op=input("please select 'c' to display coordinate sequence, 's' to display SEQRES sequence, 'a' to display the Alignment sequence, 'l' to display non water ligands in the protein or 'q' to quit:")
            if op!= 'q':
                if op.lower()in ('s','c','a','l'):
                    if op.lower()=='s':
                        """Information option-display SEQRES sequence"""
                        sef2=[]
                        rsequence2=[]
                        with open(valfile,'r')as op2:
                            for sline in op2:
                                if sline.startswith("SEQRES"):
                                    line=(str.split(sline[19:]))
                                    for rs2 in line:
                                        sef2.append(rs2)
                                        rsequence2.append(maindict[rs2])
                            seqc2=str(''.join(rsequence2).replace(' ',''))
                            Sseq = SeqRecord(Seq(seqc2,IUPAC.protein),id="RN",description=f'PDB file in {file}')
                            Sseq.format('fasta')
                            print(Sseq)
                        minioption()
                    elif op.lower()=='c':
                        """Information option-display Coordinate sequence"""
                        y=''
                        pthe=('../../python-mini-project-nanjalaruth/Results/')
                        kem1=input("please select an output file")
                        opfile1=pthe+kem1
                        with open('../../python-mini-project-nanjalaruth/Data/1zni.pdb','r')as op1:
                            with open(opfile1,'w')as q: 
                                for kline in op1:
                                    if kline.startswith("ATOM"):
                                        line=(kline[16:20])
                                        y=y+line
                                cseq = SeqRecord(Seq(str(y),IUPAC.protein),id="RN",description=f'PDB file in {file}')
                                cseq.format('fasta')
                                print(cseq)       
                        minioption()
                    elif op.lower()=='a':
                        sef=[]
                        rsequence=[]
                        with open(valfile,'r')as op:
                            for sline in op:
                                if sline.startswith("SEQRES"):
                                    line=(str.split(sline[19:70]))
                                    for rs in line:
                                        sef.append(rs)
                                        rsequence.append(maindict[rs])
                            seqc=str(''.join(rsequence).replace(' ',''))
                            seqc_list =[]
                            for i in seqc:
                                seqc_list.append(i)
                        with open(valfile,'r')as op:
                            rem = []
                            for line in op: 
                                if "REMARK 465" in line:
                                    k = line[10:15]
                                    if k == "     ":
                                        p = line[23:]
                                        if p[:4] != "    ":
                                            num= int(p)
                                            rem.append(num)
                            for q in rem:
                                seqc_list[q-1]="-X-"
                            p=("".join(seqc_list))
                            Sseq = SeqRecord(Seq(p,IUPAC.protein),id="RN",description=f'PDB file in {file}')
                            Sseq.format('fasta')
                            print(Sseq)
                        minioption()
                    elif op.lower()=='l':
                        """Information option-non water ligands in the protein"""
                        with open(valfile,'r')as op:
                            for ine in op:
                                if ine.startswith("HETATM"):
                                    if 'HOH' not in ine:
                                        print(ine)
                        minioption()
                    elif op.lower()=='q':
                        sys.exit()
                else:
                    print("Invalid option! try again")
                    minioption()
        minioption()
    else:
        print("invalid choice")
        optionI()

        
        
def optionA():
    """Alignment option"""
    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment
    
    tseq1=[]
    seq1=[]
    tseq2=[]
    seq2=[]
    
    path1= '../../python-mini-project-nanjalaruth/Results/'
    
    choice1=input("please select first PDB file you would wish to align")
    choice2=input("please select second PDB file you would wish to align")
    outputfile=input("Enter a file to write the alignment to:\n")
    
    valid_file=path+choice1
    valid_file1=path+choice2
    valid_file2=path1+outputfile
    
    if choice1 in valid_file:
        with open(valid_file,'r')as q6:
            for line in q6:
                if line.startswith("SEQRES"):
                    global a
                    a=(str.split(line[19:]))
                    for redi in a:
                        seq1.append(redi)
                        tseq1.append(main_dict[redi])
        sequence1=str(''.join(tseq1).replace(' ',''))
        Q=sequence1[:50]
        
             
    if choice2 in valid_file1:
        with open(valid_file1,'r')as q8:
            for line in q8:
                if line.startswith("SEQRES"):
                    global b
                    b=(str.split(line[19:]))
                    for rei in b:
                        seq2.append(rei)
                        tseq2.append(main_dict[rei])
        sequence2=str(''.join(tseq2).replace(' ','')) 
        R=sequence2[:50]
                        
    alignments = pairwise2.align.globalms(Q, R, 2, -1, -0.5, -0.1)
    
    with open(valid_file2,'w') as fei:
        for k in alignments:
            fei.write(format_alignment(*k))
    print(f'The alignment has successfully been written to {outputfile}')
    
    
    
    
def optionH():
    """Help option"""
    print("""
    OptionR allows you to read the contents of a file
    
    OptionS helps to look for specific sites in the PDB files that have beene read in
    
    OptionW helps to extract the content of a file and writing the output to another file
    
    OptionA helps to align sequences from two PDB files
    
    OptionI helps print out the requested information
    
    OptionH lists all option and their description in brief
    
    OptionQ terminates the program
    """)
    
    
    
def mainMenu():
    """Main menu"""
    space=" "
    global status
    status="empty"
    print('WELCOME TO:')
    print('#####'+space*3+'# #'+space*6+'######')
    print('#   #'+space*3+'#   #'+space*4+'#    #')
    print('#####'+space*3+'#    #'+space*3+'#####')
    print('#'+space*7+'#   #'+space*4+'#    #')
    print('#'+space*7+'# #'+space*6+'######')
    print(space*22+'SOFTWARE!')
    print("Protein data bank(PDB) is a repository of atomic coordinates and other information describing proteins and other important biological macromolecules")
    print("""MENU: \n
    Choices:\n
    R -Read\n 
    S- Search\n
    A- Alignment\n
    W- Write out\n
    I-Information\n
    H- Help\n
    Q- Quit\n
    """) 
    def choicesMenu():
        choice = input("Please enter 'R' to read, 'S' to search, 'A' to align, 'W' to write out, 'I' for information, 'H' for help or 'Q' to quit:")
        if choice != 'Q' and choice != "q":
            if choice.capitalize() in ('R','S','A','W','I','H'):
                if choice.capitalize() == 'R':
                    optionR()
                    choicesMenu()
                elif choice.capitalize() == 'S':
                    optionS()
                    choicesMenu()
                elif choice.capitalize() == 'A':
                    optionA()
                    choicesMenu()
                elif choice.capitalize() == 'W':
                    optionW()
                    choicesMenu()
                elif choice.capitalize() == 'I':
                    optionI()
                    choicesMenu()
                elif choice.capitalize() == 'H':
                    optionH()
                    choicesMenu()
                elif choice.capitalize() == 'Q':
                    sys.exit()
            else:
                print('invalid choice made')
                choicesMenu()
    choicesMenu()    
mainMenu()