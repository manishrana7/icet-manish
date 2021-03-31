# Ansible provisioning

## Quick start

1. [Install Ansible](https://docs.ansible.com/ansible/latest/installation_guide/intro_installation.html)
2. Prepare MacOS machine
3. Run playbook - `ansible-playbook playbook.yaml -i hosts.ini --extra-vars "ansible_sudo_pass=***"`
4. Register runner

### Prepare MacOS machine

1. Add your public SSH key to the authorized_keys
2. Install `XCode`
3. Install `brew`
4. Configure proxy

## Links

1. [Issue with building on MacOS](https://github.com/ycm-core/YouCompleteMe/issues/3770)
2. Issue with GitLab ansible module [1](https://github.com/ansible-collections/community.general/issues/440)
   , [2](https://github.com/ansible-collections/community.general/pull/1491)
3. [Issue with PATH](https://gitlab.com/gitlab-org/gitlab-runner/-/issues/3463)