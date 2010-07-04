///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

Développer un nouveau composant:
--------------------------------

1) s’inspirer de CGAL_Example dans trunk\src\components\Examples\CGAL_Example
(un renommage propre et précis tenant compte de la « case » dans les .h, .hxx et .cpp devrait vous permettre d’obtenir très facilement le squelette de votre composant)

2) le « noyau » de MEPP ne doit en principe pas être modifié

3) la description du menu propre au composant se fait dans la déclaration de la liste "actions()" (cf. fichier mepp_component_CGAL_Example_plugin.hxx comme exemple)
