# Running HippUnfold with a Vagrant VM

This option uses Vagrant to create a virtual machine that has Linux
and Singularity installed. This allows you to use Singularity to run HippUnfold 
from a clean environment, whether you are running Linux, Mac or Windows (since all
three are supported by Vagrant). Note: VirtualBox does the actual 
virtualization in this example, but Vagrant provides an easy and reproducible way to create and 
connect to the VMs (as shown below).


## Install VirtualBox and Vagrant

The example below uses Vagrant and VirtualBox installed on Ubuntu 20.04. 

The [Vagrant install instructions](https://developer.hashicorp.com/vagrant/downloads) describe
what you need to do to install on Mac, Windows or Linux.

Vagrant must use a **provider** for the actual virtualization. The instructions here assume you
are using VirtualBox for this, since it is free and easy to use, but in principle should work with any
virtualization provider.  The [VirtualBox downloads page](https://www.virtualbox.org/wiki/Downloads) 
can guide you through the process of installing it on your system (Mac, Windows, Linux supported).


## Create a Vagrant Box 

Once you have Vagrant and VirtualBox installed, the following screencast demonstrates 
how you can setup a Box with Singularity pre-loaded on it. The main steps are to 1) create a Vagrantfile, 2) start the box using `vagrant up`, and 3) connect to it using `vagrant ssh`.

Note: These screencasts are more than just videos, they are asciinema recordings -- you can pause them and then copy-paste text
directly from the asciinema cast!

```{asciinema} ../casts/vagrant_hippunfold_setup.cast
---
preload: 1
speed: 2
---
```

This is the `Vagrantfile` used in the video, for quick reference:
```
Vagrant.configure("2") do |config|
  config.vm.box = "sylabs/singularity-3.7-ubuntu-bionic64"
  config.vm.provider "virtualbox" do |vb|
     vb.cpus = 8
     vb.memory = "8096"
   end
end
```


## Download the test dataset

We are downloading the test dataset with the following 
command:
```
wget https://www.dropbox.com/s/mdbmpmmq6fi8sk0/hippunfold_test_data.tar
```


```{asciinema} ../casts/vagrant_hippunfold_get_data.cast
---
preload: 1
speed: 2
---
```


## Download the HippUnfold container

We pull/build the container from DockerHub:
```
singularity pull docker://khanlab/hippunfold:latest
```

## Run HippUnfold
This demonstrates the basic HippUnfold options, and how
to perform a dry-run:

```{asciinema} ../casts/vagrant_hippunfold_dryrun.cast
---
preload: 1
speed: 2
---
```

Finally, we can run HippUnfold using all the cores:

```{asciinema} ../casts/vagrant_hippunfold_run.cast
---
preload: 1
speed: 2
---
```


