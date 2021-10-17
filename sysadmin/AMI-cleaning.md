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
- Log in as usual - and do basic system updates
```
sudo apt update
sudo apt upgrade
sudo apt --purge autoremove
```
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
- Clean up logs, as per Pedro Lobito's comment at [https://serverfault.com/questions/185253/delete-all-of-var-log](https://serverfault.com/questions/185253/delete-all-of-var-log):
```
sudo su -
cd /var/log
set +o noclobber
for CLEAN in $(find /var/log/ -type f); do /usr/bin/cp /dev/null $CLEAN; done
```

- Remove authorized keys and history for root
```
sudo su -
shred -u .ssh/*
shred -u .*history .*hsts .gitconfig .lesshst .viminfo
rm -rf snap .cache .cpan* .vim
history -c
exit
```
- Remove authorized keys and history for ubuntu - check for other files too - careful, this will prevent any further logins to this machine
```
shred -u .ssh/*
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

Copying the public AMI to different regions (assuming you are starting in Singapore) - we're parsing with `jq` on the command line:
```
AMI_SOURCE=<AMI_ID>

SOURCE_REGION=ap-southeast-1
IMAGE_NAME=$(aws ec2 describe-images --image-ids $AMI_SOURCE | jq ".Images[].Name" | sed -e 's/^"//' -e 's/"$//')
IMAGE_DESC=$(aws ec2 describe-images --image-ids $AMI_SOURCE | jq ".Images[].Description" | sed -e 's/^"//' -e 's/"$//')
NAME_TAG=$(aws ec2 describe-images --image-ids $AMI_SOURCE | jq ".Images[].Tags[]")
# check to make sure that worked
echo $IMAGE_NAME
echo $IMAGE_DESC
echo $NAME_TAG

# 17 regions accessible by default
for REGION in us-east-1 us-east-2 us-west-1 us-west-2 ap-south-1 ap-northeast-3 ap-northeast-2 ap-southeast-2 ap-northeast-1 ca-central-1 eu-central-1 eu-west-1 eu-west-2 eu-west-3 eu-north-1 sa-east-1; do
  NEW_ID=$(aws ec2 copy-image --source-image-id $AMI_SOURCE --source-region $SOURCE_REGION --region $REGION --name "$IMAGE_NAME" | jq ".ImageId" | sed -e 's/^"//' -e 's/"$//')
  echo $REGION $NEW_ID
  aws ec2 create-tags --region $REGION --resources $NEW_ID --tags "[ $NAME_TAG ]"
done

# tagging snapshots and changing permissions needs to be done later when done
REGION=<region>
NEW_ID=<AMI ID>
NEW_SNAPSHOT=$(aws ec2 describe-images --region $REGION --image-ids $NEW_ID | jq ".Images[].BlockDeviceMappings[].Ebs.SnapshotId" | sed -e 's/^"//' -e 's/"$//')
aws ec2 create-tags --region $REGION --resources $NEW_SNAPSHOT --tags "[ $NAME_TAG ]"
aws ec2 modify-image-attribute --region $REGION --image-id $NEW_ID --launch-permission "Add=[{Group=all}]"

# Some regions require activation to access:
# af-south-1 ap-east-1 eu-south-1 me-south-1
```
