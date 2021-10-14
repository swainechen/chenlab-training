+++
title = "Local machine <-> EC2 instance"
draft = false
weight = 120
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

1.	On your local computer, open a new terminal and create a sample text file.

```bash
echo “Hello, world” > test.txt
```

2.	Use scp to copy this sample test file into the home directory of your remote machine/ instance, by running the following command on your local machine.
```bash
scp -i path/to/key.pem path/to/test.txt ubuntu@<ip-address>:/home/ubuntu/
```

If you are using windows command prompt and have the PuTTY suite of utilities installed, run the following command
```bash
pscp -P 22 -i path/to/key.ppk path/to/test.txt ubuntu@<ip-address>:/home/ubuntu/
```

Alternative: use rsync - the options here preserve timestamps and permissions (-a), and only transfer if it would result in an updated version (-u). This can help with tracking what the latest file is; using scp as above will just overwrite the destination even if the file exists at the destination and is newer.
```bash
rsync -e 'ssh -i /path/to/key.pem' -avxu /path/to/test.txt ubuntu@<ip-address>:/home/ubuntu/
```

3.	Open the test.txt file using a text editor on your remote machine and edit it - add random text and save it.

4.	Copy the test.txt file back to your local machine by using the scp/pscp command on your local terminal.

```bash
scp -i path/to/key.pem ubuntu@<ip-address>:/home/ubuntu/test.txt path/to/folder/localmachine
```
If you are using windows command prompt run the following command
```bash
pscp -P 22 -i path/to/key.ppk ubuntu@<ip-address>:/home/ubuntu/test.txt path/to/folder/localmachine
```

Alternative: use rsync:
```bash
rsync -e 'ssh -i /path/to/key.pem' -avux ubuntu@<ip-address>:/home/ubuntu/test.txt /path/to/folder/localmachine
```

5.	Open the copied file to review your edits.
