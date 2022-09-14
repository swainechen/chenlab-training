+++
title = "i. Linux commands"
date = 2021-05-23T16:24:30-04:00
draft = false 
weight = 90
tags = ["basics", "commands", "terminal"]
+++

|1. File/directory commands|Description|
|---|---|
|_**ls**_|List the files in the current folder|  
|ls -al|Same as above but the ‘a’ flag is for “all the files” and the ‘l’ flag is for long format|
|ls -trlh |List files but with ‘t’ flag to sort files by time, ‘r’ flag to sort files in reverse order, and ‘h’ flag to covert bytes to human-readable format (kB, MB, GB, etc)|   
|mkdir folder/|Create folder|   
|mkdir -p path/folder/|Create folder (creating intermediate folders as necessary, does not fail if the folder already exists)|
|cd /home/ec2-user/|Change directory to folder /home/ec2-user/|   
|mv my_file destination_folder/|Move file into destination folder|
|cp my_file ~/destination_folder/|Copy file into destination folder. The tilde ~ refers to the home directory|
|cp -ax /tmp/fastq /mnt/volume1|Same as above but with ‘a’ flag for ‘archive’ (preserve file attributes) and ‘x’ flag for ‘stay on this file system’| 
|df -h|Display free disk space. The ‘h’ flag is to convert disk size info into human-readable format (GB, TB, etc)|
|which executable_file|Shows folder location of the executable file (if it is in the path)|
|chmod 600 ~/.ssh/id_rsa|Set the ssh key permissions so that only the file owner can read and write it|
|chmod 600 my_ssh_key.pem|Same as above but applied to my_ssh_key.pem|   
|||   
|**2. File Transfer Commands**|
|wget https://hostname/path/file_to_download|Download file at the URL specified|
|curl -LO https://hostname/path/file_to_download|Download file at the URL specified ‘L’ flag to handle redirects, ‘O’ flag to use the last bit of the URL after the slash as the filename| 	
|scp my_file ubuntu@32.18.231.125:/home/ubuntu/|Copy the file my_file to the server at 32.18.231.125, at folder /home/ubuntu/, logging in as the user ‘ubuntu’ (you will be prompted for password)|
|scp -I my_ssh_keypair.pem my_file ubuntu@32.18.231.125:/home/ubuntu/|Same as above but specifying the keypair as my_ssh_keypair.pem|
|pscp -P 22 -I my_ssh_keypair.ppk my_file ubuntu@32.18.231.125:/home/ubuntu/|Same as above but using pscp (PuTTY secure copy client)|
|rsync -e 'ssh -I my_ssh_keypair.pem’ -avxu my_file ubuntu@32.18.231.125:/home/ubuntu/|Same as above but using rsync|
|Shortcuts (requires configuring ~/.ssh/config)||
|scp my_file myserver:/home/ubuntu/||
|rsync -avxu my_file my_server:/home/ubuntu/||
|Configuration of ~/.ssh/config with hostname, username and key pair:  	
"Host myserver
    Hostname 32.18.231.125
    User ubuntu
    IdentityFile ~/.ssh/my_ssh_keypair.pem"||   
|||
|**3. SSH (Secure Shell) Commands**|	
|ssh ec2-user@32.18.231.125|Uses SSH to log in to the server at IP address 32.18.231.125, with the username “ec2-user”|
|ssh -i my_ssh_keypair.pem ec2-user@32.18.231.125|Same as above but specifying my_ssh_keypair.pem as the keypair (for password-less login)|
|ssh ec2-user@32.18.231.125 “ls -al”|Same as first line, but runs command “ls -al” on the server and then logs out|  
|||     	
|**4. Mounting a volume, formatting it, and then unmounting it**|	
|lsblk|List block devices|
|sudo file -s /dev/nvme1n1|Check if the volume has any data|
|sudo mkfs -t ext4 /dev/nvme1n1|Format the volume to the ext4 filesystem|
|sudo mkdir /mnt/volume1|Create a folder where the volume will be mounted|
|sudo mount /dev/nvme1n1 /mnt/volume1|Mount the volume “/dev/nvme1n1” to the folder “/mnt/volume1”|
|sudo chmod -R 777 /mnt/volume1|Set the permissions of the file volume at /mnt/volume1 to 777. The ‘R’ flag is to recursively do this for all subfolders|
|sudo chown -R ubuntu /mnt/volume1|Change ownership of the drive so that the user “ubuntu” can access it (and not just the root user). The ‘R’ flag is to recursively do this for all subfolders|
|sudo umount /dev/nvme1n1|Unmount the volume|    
|||    	
|**5. Software Package Installation**|
|***Amazon Linux, CentOS, & Red Hat***|
|sudo yum install nmap|Install nmap|
|yum install gcc git zlib-devel|Install the packages “gcc”, “git”, and “zlib-devel” (will not work if you are not the root user)|
|***Ubuntu & Debian***|
|sudo apt update && sudo apt -y upgrade|Update the system, then upgrade|
|sudo apt install nmap|Install nmap|
|dpkg –list \| grep awscli|List all the currently installed packages, the pipe it through "grep" and print only those whose name contains "awscli"|
||| 	
|***Python***|	
|pip3 list\| grep awscli|Check whether the awscli Python package is installed|
||| 	
|***Other commands for software package installation***|
|git clone https://github.com/lh3/seqtk.git|Clone the Github repo “seqtk” into a directory “seqtk/” inside the current folder|
|gunzip my_program.zip|Unzip the contents of .zip file into the current folder|
|make|A command to compile code (requires a Makefile in the same folder for instructions)|
||| 	
|***Root privilege escalation and switching accounts***|
|sudo su -|Switch to root user account (the "-" option makes sure your environment and home directory also switch to the root account's versions)|
|sudo su username|Switch to username account|
|sudo [command]|Run command as the root user|    
|||     	
|**6. Genomic-specific command-line-based programs**|	
|bwa|http://bio-bwa.sourceforge.net/|
|lofreq|https://csb5.github.io/lofreq/|
|minimap2|https://github.com/lh3/minimap2|
|samtools|https://github.com/samtools/samtools|    
|||    	
|**7. AWS CLI**|
|aws help|Get help on the AWS CLI|
|aws configure|Set up credentials for the AWS CLI (requires AWS Access Key ID and AWS Secret Access Key)|
|aws ec2 help|Get help on the AWS CLI / EC2 command|
|aws ec2 describe-instances|Describe EC2 instances (returns a JSON)|
|aws ec2 describe-key-pairs|Describe EC2 key pairs (returns a JSON)|
|aws s3 mb s3://my_bucket_name|Make S3 bucket with the name my_bucket_name|
|aws s3 ls s3://my_bucket_name|Lists the contents of the S3 bucket|
|aws s3 ls –profile myprofilename s3://my_bucket_name|Same as above but uses profile “myprofilename”|
|aws s3 cp help|Get help on the AWS CLI / S3 / CP command|
|aws s3 cp s3://mybucket/myobject ./|Copies myobject in the S3 bucket “mybucket” to current folder|
|aws s3 cp ./my_file_to_upload s3://my_bucket_name|Copies my_file_to_upload to S3 bucket|
|aws s3 cp novel_annotations.UCSC.gtf s3://mybucket/path/ --acl public-read|Copies file to S3 bucket and sets the ACL to public read|
|aws s3 sync –no-sign-request s3://sg-nex-data/data/bambu_training/bam/ ./|Sync S3 folder to current directory|   
