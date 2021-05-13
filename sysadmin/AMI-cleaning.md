# Prepping an AMI for release
Some cleaning needs to be done prior to releasing a public AMI.
To get it on AWS Marketplace, it also needs to be in us-east-1 (?)

Some resources:
* https://docs.aws.amazon.com/marketplace/latest/userguide/best-practices-for-building-your-amis.html
* https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/building-shared-amis.html
* https://docs.aws.amazon.com/marketplace/latest/userguide/product-and-ami-policies.html

Brief steps:
- First make an initial AMI - I call this pre-public
- Then launch an instance off this pre-public AMI
- Log in as usual
- Disable password-based remote logins for root
- Edit /etc/ssh/sshd_config, need this line:
```
PermitRootLogin without-password
```
Should also check that no password ssh logins are allowed:
```
PasswordAuthentication no
```
- Disable local root access
```
sudo passwd -l root
```
- Remove SSH host key pairs
```
sudo shred -u /etc/ssh/*_key /etc/ssh/*_key.pub
```
- Remove authorized keys - careful, this will prevent any further logins to this machine
```
sudo su -
shred -u .ssh/*
exit
# now as ubuntu
shred -u .ssh/*
```
- Remove history - check for other files too
```
sudo su -
shred -u .*history .*hsts .gitconfig .lesshst .viminfo
rm -rf snap .cache .cpan* .vim
history -c
exit
# now as ubuntu
shred -u .*history .*hsts .gitconfig .lesshst .viminfo .sudo_as_admin_successful
rm -rf .cache .config .local .mozilla .vim .keras
history -c
exit
```
- Now stop the instance on the console, then make an image.
- Be sure to label the AMI itself and the snapshot backing it.
- Relaunch and test logging in, software, etc.
- Try the self-service scan: https://aws.amazon.com/marketplace/management/manage-products
- Then release

