opemidimejiosatoyinbo@cloudshell:~/biocomputing$ # Going back to the root folder to activate conda, then I open a new terminal to launch
opemidimejiosatoyinbo@cloudshell:~/biocomputing$ cd ../
opemidimejiosatoyinbo@cloudshell:~$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

opemidimejiosatoyinbo@cloudshell:~$ # Activating conda environment
opemidimejiosatoyinbo@cloudshell:~$ bash Miniconda3-latest-Linux-x86_64.sh

# In order to continue the installation process, I reviewed the license agreement (all od these were done in CloudShell before opening a new terminal to install Bioinformatics softwares).

(base) opemidimejiosatoyinbo@cloudshell:~$ # Installing Bioinformatics Software on the Terminal
(base) opemidimejiosatoyinbo@cloudshell:~$ #--------------------
(base) opemidimejiosatoyinbo@cloudshell:~$ # checking conda version
(base) opemidimejiosatoyinbo@cloudshell:~$ conda --version
conda 25.7.0
(base) opemidimejiosatoyinbo@cloudshell:~$ #--------------------
(base) opemidimejiosatoyinbo@cloudshell:~$ # 1. Activate your base conda environment
(base) opemidimejiosatoyinbo@cloudshell:~$ conda activate base
(base) opemidimejiosatoyinbo@cloudshell:~$ #--------------------
(base) opemidimejiosatoyinbo@cloudshell:~$ # Configuring bioconda channels
(base) opemidimejiosatoyinbo@cloudshell:~$ conda config --add channels defaults
(base) opemidimejiosatoyinbo@cloudshell:~$ conda config --add channels bioconda
(base) opemidimejiosatoyinbo@cloudshell:~$ conda config --add channels conda-forge
(base) opemidimejiosatoyinbo@cloudshell:~$ conda config --set channel_priority strict
(base) opemidimejiosatoyinbo@cloudshell:~$ #--------------------
(base) opemidimejiosatoyinbo@cloudshell:~$ # 2. Create a conda environment named 'funtools'
(base) opemidimejiosatoyinbo@cloudshell:~$ conda create -n funtools python=3.10 -y                                                                   
(base) opemidimejiosatoyinbo@cloudshell:~$ #--------------------                                         
(base) opemidimejiosatoyinbo@cloudshell:~$ # 3. Activate the 'funtools' environment
(base) opemidimejiosatoyinbo@cloudshell:~$ conda activate funtools                                       
(funtools) opemidimejiosatoyinbo@cloudshell:~$ #--------------------                                     
(funtools) opemidimejiosatoyinbo@cloudshell:~$ # 4. Install Figlet using conda or apt-get
(funtools) opemidimejiosatoyinbo@cloudshell:~$ # Using conda:
(funtools) opemidimejiosatoyinbo@cloudshell:~$ conda install -y -c conda-forge figlet                    
(funtools) opemidimejiosatoyinbo@cloudshell:~$ # Package was not found using conda, so I will proceed to 
use apt-get: 
(funtools) opemidimejiosatoyinbo@cloudshell:~$ sudo apt-get install -y figlet
(funtools) opemidimejiosatoyinbo@cloudshell:~$ #--------------------
(funtools) opemidimejiosatoyinbo@cloudshell:~$ # 5. Run 'figlet' <my name>
(funtools) opemidimejiosatoyinbo@cloudshell:~$ figlet "Opemidimeji"
  ___                       _     _ _                 _ _ 
 / _ \ _ __   ___ _ __ ___ (_) __| (_)_ __ ___   ___ (_|_)
| | | | '_ \ / _ \ '_ ` _ \| |/ _` | | '_ ` _ \ / _ \| | |
| |_| | |_) |  __/ | | | | | | (_| | | | | | | |  __/| | |
 \___/| .__/ \___|_| |_| |_|_|\__,_|_|_| |_| |_|\___|/ |_|
      |_|                                          |__/   
(funtools) opemidimejiosatoyinbo@cloudshell:~$ #--------------------
(funtools) opemidimejiosatoyinbo@cloudshell:~$ # 6. Install 'bwa' through the bioconda channel
(funtools) opemidimejiosatoyinbo@cloudshell:~$ conda install -y -c bioconda bwa
(funtools) opemidimejiosatoyinbo@cloudshell:~$ #--------------------
(funtools) opemidimejiosatoyinbo@cloudshell:~$ # 7. Install blast through the bioconda channel
(funtools) opemidimejiosatoyinbo@cloudshell:~$ conda install -y -c bioconda blast                                                                      
(funtools) opemidimejiosatoyinbo@cloudshell:~$ #--------------------                                     
(funtools) opemidimejiosatoyinbo@cloudshell:~$ # 8. Install 'samtools' through the bioconda channel
(funtools) opemidimejiosatoyinbo@cloudshell:~$ conda install -y -c bioconda samtools                                                                                                  
(funtools) opemidimejiosatoyinbo@cloudshell:~$ #--------------------
(funtools) opemidimejiosatoyinbo@cloudshell:~$ # 9. Install 'bedtools' through the bioconda channel
(funtools) opemidimejiosatoyinbo@cloudshell:~$ conda install -y -c bioconda bedtools
(funtools) opemidimejiosatoyinbo@cloudshell:~$ #--------------------
(funtools) opemidimejiosatoyinbo@cloudshell:~$ # 10. Install 'spades.py' through the bioconda channel
(funtools) opemidimejiosatoyinbo@cloudshell:~$ conda install -y -c bioconda spades
(funtools) opemidimejiosatoyinbo@cloudshell:~$ #--------------------
(funtools) opemidimejiosatoyinbo@cloudshell:~$ # 11. Install 'bcftools' through the bioconda channel
(funtools) opemidimejiosatoyinbo@cloudshell:~$ conda install -y -c bioconda bcftools
(funtools) opemidimejiosatoyinbo@cloudshell:~$ #--------------------
(funtools) opemidimejiosatoyinbo@cloudshell:~$ # 12. Install 'fastp' through the bioconda channel
(funtools) opemidimejiosatoyinbo@cloudshell:~$ conda install -y -c bioconda fastp                                                                             
(funtools) opemidimejiosatoyinbo@cloudshell:~$ #--------------------                                     
(funtools) opemidimejiosatoyinbo@cloudshell:~$ # 13. Install 'multiqc' through the bioconda channel      
(funtools) opemidimejiosatoyinbo@cloudshell:~$ conda install -y -c bioconda multiqc                      

(funtools) opemidimejiosatoyinbo@cloudshell:~$ 