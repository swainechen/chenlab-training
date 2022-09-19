+++
title = "i. Linux commands"
date = 2021-05-23T16:24:30-04:00
draft = false 
weight = 90
tags = ["basics", "commands", "terminal"]
+++

|1. FILE/DIRECTORY COMMANDS|DESCRIPTION|
|---|---|
|``` ls ```|_List the files in the current folder_|
|```ls -al```|_Same as above but the ‘a’ flag is for “all the files” and the ‘l’ flag is for long format_|
|```ls -trlh```|_List files but with ‘t’ flag to sort files by time, ‘r’ flag to sort files in reverse order, and ‘h’ flag to covert bytes to human-readable format (kB, MB, GB, etc)_|   
|```mkdir folder/```|_Create folder_|   
|```mkdir -p path/folder/```|_Create folder (creating intermediate folders as necessary, does not fail if the folder already exists)_|
|```cd /home/ec2-user/```|_Change directory to folder /home/ec2-user/_|   
|```mv my_file destination_folder/```|_Move file into destination folder_|
|```cp my_file ~/destination_folder/```|_Copy file into destination folder. The tilde ~ refers to the home directory_|
|```cp -ax /tmp/fastq /mnt/volume1```|_Same as above but with ‘a’ flag for ‘archive’ (preserve file attributes) and ‘x’ flag for ‘stay on this file system’_| 
|```df -h```|_Display free disk space. The ‘h’ flag is to convert disk size info into human-readable format (GB, TB, etc)_|
|```which executable_file```|_Shows folder location of the executable file (if it is in the path)_|
|```chmod 600 ~/.ssh/id_rsa```|_Set the ssh key permissions so that only the file owner can read and write it_|
|```chmod 600 my_ssh_key.pem```|_Same as above but applied to my_ssh_key.pem_|   
|||  
||| 
|**2. FILE TRANSFER COMMANDS**|
|```wget https://hostname/path/file_to_download```|_Download file at the URL specified_|
|```curl -LO https://hostname/path/file_to_download```|_Download file at the URL specified ‘L’ flag to handle redirects, ‘O’ flag to use the last bit of the URL after the slash as the filename_| 	
|```scp my_file ubuntu\@32.18.231.125:/home/ubuntu/```|_Copy the file my_file to the server at 32.18.231.125, at folder /home/ubuntu/, logging in as the user ‘ubuntu’ (you will be prompted for password)_|
|```scp -I my_ssh_keypair.pem my_file ubuntu\@32.18.231.125:/home/ubuntu/```|_Same as above but specifying the keypair as my_ssh_keypair.pem_|
|```pscp -P 22 -I my_ssh_keypair.ppk my_file ubuntu\@32.18.231.125:/home/ubuntu/```|_Same as above but using pscp (PuTTY secure copy client)_|
|```rsync -e 'ssh -I my_ssh_keypair.pem’ -avxu my_file ubuntu\@32.18.231.125:/home/ubuntu/```|_Same as above but using rsync_|
|**~ Shortcuts (requires configuring ~/.ssh/config) ~**||
|```scp my_file myserver:/home/ubuntu/```||
|```rsync -avxu my_file my_server:/home/ubuntu/```||
|**~ Configuration of ~/.ssh/config with hostname, username and key pair: ~**  	
Host: _myserver_
Hostname: _32.18.231.125_
User: _ubuntu_
IdentityFile: _~/.ssh/my_ssh_keypair.pem_|   
|||
|||
|**3. SSH (SECURE SHELL) COMMANDS**|	
|```ssh ec2-user\@32.18.231.125```|_Uses SSH to log in to the server at IP address 32.18.231.125, with the username “ec2-user”_|
|```ssh -i my_ssh_keypair.pem ec2-user\@32.18.231.125```|_Same as above but specifying my_ssh_keypair.pem as the keypair (for password-less login)_|
|```ssh ec2-user\@32.18.231.125 “ls -al”```|_Same as first line, but runs command “ls -al” on the server and then logs out_|  
||| 
||||    	
|**4. MOUNT, FORMAT, UNMOUNT A VOLUME**|	
|```lsblk```|_List block devices_|
|```sudo file -s /dev/nvme1n1```|_Check if the volume has any data_|
|```sudo mkfs -t ext4 /dev/nvme1n1```|_Format the volume to the ext4 filesystem_|
|```sudo mkdir /mnt/volume1```|_Create a folder where the volume will be mounted_|
|```sudo mount /dev/nvme1n1 /mnt/volume1```|_Mount the volume “/dev/nvme1n1” to the folder “/mnt/volume1”_|
|```sudo chmod -R 777 /mnt/volume1```|_Set the permissions of the file volume at /mnt/volume1 to 777. The ‘R’ flag is to recursively do this for all subfolders_|
|```sudo chown -R ubuntu /mnt/volume1```|_Change ownership of the drive so that the user “ubuntu” can access it (and not just the root user). The ‘R’ flag is to recursively do this for all subfolders_|
|```sudo umount /dev/nvme1n1```|_Unmount the volume_|    
||| 
|||   	
|**5. SOFTWARE PACKAGE INSTALLATION**|
|~ ***Amazon Linux, CentOS, & Red Hat*** ~|
|```sudo yum install nmap```|_Install nmap_|
|```yum install gcc git zlib-devel```|_Install the packages “gcc”, “git”, and “zlib-devel” (will not work if you are not the root user)_|
|~ ***Ubuntu & Debian*** ~|
|```sudo apt update && sudo apt -y upgrade```|_Update the system, then upgrade_|
|```sudo apt install nmap```|_Install nmap_|
|```dpkg –list \| grep awscli```|_List all the currently installed packages, the pipe it through "grep" and print only those whose name contains "awscli"_|
||| 	
|~ ***Python*** ~|	
|```pip3 list\| grep awscli```|_Check whether the awscli Python package is installed_|
||| 	
|~ ***Other commands for software package installation*** ~|
|```git clone https://github.com/lh3/seqtk.git```|_Clone the Github repo “seqtk” into a directory “seqtk/” inside the current folder_|
|```gunzip my_program.zip```|_Unzip the contents of .zip file into the current folder_|
|```make```|_A command to compile code (requires a Makefile in the same folder for instructions)_|
||| 	
|~ ***Root privilege escalation and switching accounts*** ~|
|```sudo su -```|_Switch to root user account (the "-" option makes sure your environment and home directory also switch to the root account's versions)_|
|```sudo su username```|_Switch to username account_|
|```sudo [command]```|_Run command as the root user_|    
||| 
|||    	
|**6. GENOME-SPECIFIC COMMAND-LINE BASED PROGRAMS**|	
|bwa|_http://bio-bwa.sourceforge.net/_|
|lofreq|_https://csb5.github.io/lofreq/_|
|minimap2|_https://github.com/lh3/minimap2_|
|samtools|_https://github.com/samtools/samtools_|    
||| 
|||   	
|**7. AWS CLI**|
|```aws help```|_Get help on the AWS CLI_|
|```aws configure```|_Set up credentials for the AWS CLI (requires AWS Access Key ID and AWS Secret Access Key)_|
|```aws ec2 help```|_Get help on the AWS CLI / EC2 command_|
|```aws ec2 describe-instances```|_Describe EC2 instances (returns a JSON)_|
|```aws ec2 describe-key-pairs```|_Describe EC2 key pairs (returns a JSON)_|
|```aws s3 mb s3://my_bucket_name```|_Make S3 bucket with the name my_bucket_name_|
|```aws s3 ls s3://my_bucket_name```|_Lists the contents of the S3 bucket_|
|```aws s3 ls –profile myprofilename s3://my_bucket_name```|_Same as above but uses profile “myprofilename”_|
|```aws s3 cp help```|_Get help on the AWS CLI / S3 / CP command_|
|```aws s3 cp s3://mybucket/myobject ./```|_Copies myobject in the S3 bucket “mybucket” to current folder_|
|```aws s3 cp ./my_file_to_upload s3://my_bucket_name```|_Copies my_file_to_upload to S3 bucket_|
|```aws s3 cp novel_annotations.UCSC.gtf s3://mybucket/path/ --acl public-read```|_Copies file to S3 bucket and sets the ACL to public read_|
|```aws s3 sync –no-sign-request s3://sg-nex-data/data/bambu_training/bam/ ./```|_Sync S3 folder to current directory_|   
|||
|||
|**8. OTHER REFERENCES**|
https://ubuntu.com/tutorials/command-line-for-beginners
https://www.digitalocean.com/community/tutorial_series/getting-started-with-linux
https://www.freecodecamp.org/news/the-linux-commands-handbook
