from Bio.Seq import Seq
import threading
import PySimpleGUI as sg
import pandas as pd
import timeit
from Bio.Blast.Applications import NcbiblastpCommandline 
import concurrent.futures
import regex as re

def logic_thread(values, window):
    
    start = timeit.default_timer()
    DNA = Seq("".join(l.rstrip() for l in values.get("-IN-")))

    def transl(DNA):
        frame1= DNA[0:].translate() 
        frame2= DNA[1:].translate()
        frame3= DNA[2:].translate()
        with concurrent.futures.ThreadPoolExecutor () as executor:
            results = [executor.submit(findORFs, frame1, 1),executor.submit(findORFs, frame2, 2),executor.submit(findORFs, frame3, 3)]
            try:
                f= open ("all_ORFs.csv", "x")
                for result in results:
                    for key, value in result.result().items():
                            f.write(">"+key + "\n")
                            f.write(value + "\n")
                f.close()
            except:
                f= open ("all_ORFs.csv", "w")
                for result in results:
                    for key, value in result.result().items():
                            f.write(">"+key + "\n")
                            f.write(value + "\n")
                f.close()
        window["-LOG1-"].update("Selecting ORFs **DONE**, BLASTing now, this might take a while...", text_color="White")
        BLAST()

    

    def findORFs(frame, framenumber):
        ORFs = {}
        hitcounter = 1
        hits = re.finditer(r"M.*?(?=\*|$)", str(frame), flags=re.DOTALL)
        for match in hits:
            s = match.start()
            e = match.end()
            if e-s >= 40:
                DNAstart= 3*s+framenumber
                DNAend= 3*e+framenumber+2
                complexcounter = str(framenumber)+"."+str(hitcounter)+"_"+str(DNAstart)+"-"+str(DNAend)
                hitcounter +=1
                ORFs[complexcounter] = match.group(0)
        return ORFs

    def BLAST():

        cline = NcbiblastpCommandline(cmd="blastp", query = "all_ORFs.csv" , db = "sum_phage_data_db", outfmt ="6", out = "results.csv") # Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
        cline()
        
        blastdf = pd.read_csv("results.csv", delimiter="\t", low_memory=False, header=None)
        blastdf.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend","sstart", "send", "evalue", "bitscore"]
        blastdf.replace(to_replace=r"'", value=" ", regex=True, inplace=True)

        selecteddf = blastdf.groupby("qseqid").head(1)
        selecteddf.to_csv("best_hits_selected.csv", sep="\t",header=True, index=False )
        
        stop = timeit.default_timer()

        window.write_event_value('-LOG2-', str(round(stop-start)))
        window["-OUTPUT-"].update(selecteddf[["qseqid","sseqid","bitscore"]].to_string())

    transl(DNA)
    

def the_gui():

    sg.theme('DarkGreen4')

    layout = [
                
                [sg.Text('Please insert your sequence here:', key= '-LOG2-')],
                [sg.Multiline(size=(50,10), key='-IN-',expand_x=(500),expand_y=(500),rstrip=True)],
                [sg.Button('Run ORF Finder', button_color='Grey',bind_return_key=True), sg.Button('Reset', button_color='Grey'), sg.Button('Exit', button_color='Grey')],
                [sg.Text(key='-LOG1-')],
                [sg.Multiline('This is your OUTPUT window.', size=(50,10), key='-OUTPUT-',expand_x=(500),expand_y=(500),rstrip=True, auto_size_text=True,)],
          ]

    window=sg.Window('ORF Finder',layout,resizable=True,finalize=True, auto_size_text=True,grab_anywhere=True)


    while True:
        event, values = window.read()
        if event in (sg.WIN_CLOSED, 'Exit'):
            break
        elif event == 'Run ORF Finder' and not values['-IN-']:
            window['-LOG1-'].update("MUST INSERT SEQUENCE TO PROCEED!", text_color="Red")
        elif event == 'Run ORF Finder' :
            window['-OUTPUT-']('')
            threading.Thread(target=logic_thread, args=(values,window), daemon=True).start()
        elif event == 'Reset':
            window['-IN-']('')
            window['-LOG1-']('')
            window['-OUTPUT-']('This is your OUTPUT window.')
        elif event == '-LOG2-':
            window['-LOG1-']("ORF Finder ** DONE ** in: " + " " + values[event] + " s")

    window.close()

if __name__ == '__main__':
    the_gui()
