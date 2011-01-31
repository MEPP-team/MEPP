///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010-2011
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

Sous Windows :
--------------

1) installation et configuration de TortoiseSVN (Sub)Version Control:

- télécharger la dernière version en date pour Windows à l’adresse suivante : http://tortoisesvn.net/downloads.html
- installer le logiciel et redémarrer tout de suite (important)

- pour démarrer facilement avec TortoiseSVN, je vous conseille de lire ce petit tutoriel: http://www.jonathan-petitcolas.com/fr/utiliser-subversion-svn-avec-tortoise
- pour une utilisation plus approfondie de TortoiseSVN, je vous recommande ce document: http://kevin.fardel.perso.esil.univmed.fr/documentation/TutorielTortoiseSVN.pdf

2) récupération du projet MEPP:

- créer un dossier sur votre système qui hébergera votre version de MEPP (ex. : MEPP_SVN) : c’est dans ce dossier que nous allons procéder au "checkout" du projet
- faire un "bouton-droit" sur le dossier créé ci-dessus puis choisir "SVN Checkout" dans le menu contextuel
- dans la fenêtre qui s’ouvre, renseigner le champ "URL of repository" par l’url du projet MEPP sur la forge : https://nom-du-développeur@scm.gforge.liris.cnrs.fr/svnroot/mepp
en prenant bien soin de renseigner "votre_username_gforge" puis valider par "Ok"
- à la première utilisation vous pouvez obtenir une alerte de sécurité relatif à la clef SSH (et/ou à un certificat), répondre "Oui" pour stocker la clef dans le cache
- une demande d'authentification vous est ensuite demandée
- le "checkout" commence alors

3) abonnement aux listes de diffusion: https://gforge.liris.cnrs.fr/mail/?group_id=39


Sous Linux :
------------

1) installation et configuration de SVN (Sub)Version Control:

- utiliser la commande "sudo apt-get install subversion" ou le "gestionnaire de paquets synaptic" pour installer le paquet "subversion"

2) récupération du projet MEPP:

- créer un dossier sur votre système qui hébergera votre version de MEPP (ex. : MEPP_SVN) : c’est dans ce dossier que nous allons procéder au "checkout" du projet
- utiliser la commande "svn checkout https://nom-du-développeur@scm.gforge.liris.cnrs.fr/svnroot/mepp" en prenant bien soin de renseigner "votre_username_gforge"
- à la première utilisation vous pouvez obtenir une alerte de sécurité relatif à la clef SSH (et/ou à un certificat), répondre "yes" pour stocker la clef dans le cache
- une demande d'authentification vous est ensuite demandée
- le "checkout" commence alors

3) abonnement aux listes de diffusion: https://gforge.liris.cnrs.fr/mail/?group_id=39


Sous Mac OS X :
---------------

1) installation et configuration de SVN (Sub)Version Control:

- installer Xcode: http://developer.apple.com/technologies/xcode.html
- installer MacPorts: http://www.macports.org/install.php
- mettre à jour MacPorts: sudo port -v selfupdate
- installer le paquet "subversion" avec MacPorts: sudo port install subversion

2) récupération du projet MEPP:

- créer un dossier sur votre système qui hébergera votre version de MEPP (ex. : MEPP_SVN) : c’est dans ce dossier que nous allons procéder au "checkout" du projet
- utiliser la commande "svn checkout https://nom-du-développeur@scm.gforge.liris.cnrs.fr/svnroot/mepp" en prenant bien soin de renseigner "votre_username_gforge"
- à la première utilisation vous pouvez obtenir une alerte de sécurité relatif à la clef SSH (et/ou à un certificat), répondre "yes" pour stocker la clef dans le cache
- une demande d'authentification vous est ensuite demandée
- le "checkout" commence alors

3) abonnement aux listes de diffusion: https://gforge.liris.cnrs.fr/mail/?group_id=39
