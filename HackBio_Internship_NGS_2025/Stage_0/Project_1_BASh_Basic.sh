  197  #!/bin/bash
  198  #--------------------
  199  # Biocomputing Task
  200  # By: Opemidimeji Osatoyinbo
  201  # Date: 30/08/2025
  202  # BASh Basic and Installing Bioinformatics Software on the Terminal
  203  #--------------------
  204  # 1. Print your name
  205  echo "Opemidimeji"
  206  #--------------------
  207  # 2. Create a folder titled your name
  208  mkdir Opemidimeji
  209  #--------------------
  210  # 3. Create another new directory titled 'biocomputing' and change to that directory with one line command
  211  mkdir biocomputing
  212  cd biocomputing
  213  #--------------------
  214  # 4. Downloading files
  215  wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
  216  wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
  217  wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk.1
  218  #--------------------
  219  # 5. Moving the  '.fna' file to the folder titled my name
  220  mv wildtype.fna ../Opemidimeji/
  221  #--------------------
  222  # 6. Delete the duplicate '.gbk' file
  223  rm wildtype.gbk.1
  224  #--------------------
  225  # 7. Confirm if the '.fna' file is mutant or wild type ('tatatata' vs 'tata')
  226  # Mutant = contains repeat expansion 'tatatata'
  227  # Wildtype = contains shorter 'tata' repeat
  228  # change directory back to ../Opemidimeji/ as it contains the '.fna' file
  229  cd ../Opemidimeji
  230  check_mutation () { local wildype.fna=$1; if grep "tatatata" "$wildtype.fna"; then echo "Mutant sequence detected in $wildtype.fna"; grep "tatatata" "$wildtype.fna" > mutant_sequence.txt; echo "Mutant lines save to mutant_sequence.txt"; else echo "Wildtype sequence in $wildtype.fna"; fi; }
  231  grep "tatatata" wildtype.fna
  232  #--------------------
  233  # 8. if mutant, print all matching lines into a new file
  234  grep "tatatata" ../Opemidimeji/wildtype.fna > mutant_sequence.txt
  235  #--------------------
  236  # 9. Count number of lines (excluding header) in the '.gbk' file
  237  # change directory to biocomputing containing the '.gbk' file
  238  cd ../biocomputing
  239  grep -v "^LOCUS" wildtype.gbk | wc -l
  240  #--------------------
  241  # 10. Print the sequence length of the '.gbk' file. (Use the 'LOCUS' tag in the firt line
  242  grep "^LOCUS" wildtype.gbk | awk '{print $3, $4}'
  243  #--------------------
  244  # 11. Print the source organism of the '.gbk' file. (Use the SOURCE' tag in the first line)
  245  grep "SOURCE" wildtype.gbk
  246  #--------------------
  247  # 12. List all the gene names of the '.gbk' file. Hint {'grep '/gene=''}
  248  grep "/gene=" wildtype.gbk
  249  # 13. Clear your terminal space and print all command used today
  250  clear
  251  history
opemidimejiosatoyinbo@cloudshell:~$ #---------------------
opemidimejiosatoyinbo@cloudshell:~$ # 14. List the files in the two folders
opemidimejiosatoyinbo@cloudshell:~/Opemidimeji$ #I had to open the directory 'Opemidimeji' before listing
opemidimejiosatoyinbo@cloudshell:~$ cd ../Opemidimeji
opemidimejiosatoyinbo@cloudshell:~/Opemidimeji$ ls
mutant_sequence.txt  wildtype.fna
opemidimejiosatoyinbo@cloudshell:~/Opemidimeji$ #----------------------
opemidimejiosatoyinbo@cloudshell:~/Opemidimeji$ # Moving to 'biocomputing' before listing too
opemidimejiosatoyinbo@cloudshell:~/Opemidimeji$ cd ../biocomputing
opemidimejiosatoyinbo@cloudshell:~/biocomputing$ ls
wildtype.gbk