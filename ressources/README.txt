# Explication sur le script job_clingo.sh
La première ligne indique que le script doit être lu avec l'interpréteur bash.
Les entêtes #SBATCH indiquent que le reste de la ligne est un argument de la commande sbatch. Voici le détail des options présentes dans le script:

## Options de la commande sbatch

-t : le temps alloué au travail. J'ai l'habitude d'écrire le nombre d'heures total (format `HH:MM:SS`), mais sachez qu'inti accepte une écriture de type `N-HH:MM:SS` pour "N jours". Sans cette option, le temps par défaut est de 1h.
/!\ Le temps de calcul maximum possible dépend de la combinaison partition/queue de travail: cf options `-p` et `--qos`.

-N : nombre de noeuds (nodes). Un noeud, c'est une machine. Vous n'en aurez besoin que d'une seule, a priori.

-n : nombre de tâches à exécuter en parallèle. Encore une fois, ici c'est 1 : la tâche en elle-même est parallélisée, mais ça reste une unique tâche.

-c : nombre de coeurs (cpus) utilisés. Je l'ai mis à 8, mais vous pouvez l'augmenter.

--mem : mémoire demandée pour votre travail. Si le logiciel consomme plus de mémoire que ce que vous en avez demandé, l'ordonnanceur va faire planter votre travail. Vous pouvez la modifier. Cette option accepte des dénominations comme "M" (pour Mo, par défaut) ou G (pour Go).

-p : le nom de la partition. Vous pouvez demander normal, xlarge ou xxlarge (cf temps de calcul et qos, plus bas).

--qos : queue de travail. Il y en a plusieurs susceptibles de vous intéresser:
1) "default" sur la partition "normal" vous permet de lancer des travaux jusqu'à 24h;
2) "long" sur la parition "normal" vous permet de lancer des travaux jusqu'à 72h;
3) "xlarge"/"xxlarge" sur les partitions "xalrge" ou "xxlarge" respectivement, vous permet de lancer des travaux pour 7 jours;
4) "xlarge_month"/"xxlarge_month" vous permet de lancer des travaux pour 31 jours.

Les qos 1) et 2) sont exclusives à la partition "normal", et les qos 3) et 4) sont exclusives à "xlarge"/"xxlarge".

--job-name : nom du job indiqué sur le flux de travail. Mettez ce que vous voulez.

--mail-user : l'adresse mail où envoyer les mails de changement d'état (cf ci-dessous)
--mail-type : envoie un mail lorsque le travail change d'état dans le flux de travail. Les possibilités sont:
BEGIN : quand le travail est envoyé sur un noeud et débute réellement ;
FAIL : quand le travail plante ;
END : quand le travail termine ;
ALL : les trois options précédentes.

L'envoi de mail est facultatif. Si vous ne désirez pas recevoir de mail d'inti, effacez ces deux lignes.
Ecrivez bien votre adresse mail, au risque de vous faire remonter les bretelles (je parle d'expérience personnelle).
 
-D : chemin du répertoire de travail. Mettez ce que vous voulez.
/!\ TOUS les chemins relatifs que vous mettrez dans ce script seront définis à partir du répertoire de travail. En cas de doute, écrivez un chemin absolu.

--error: fichier dans lequel sera inscrit l'erreur standard du travail. Le format dans le script que je vous envoie est de type <nom_du_job>_<id_job>-<id_utilisateur>.err

--output: fichier où sera écrite la sortie standard du travail. Le format est de type <nom_du_job>_<id_job>-<id_utilisateur>.out

Vous pouvez combiner erreur et sortie dans un même fichier, surtout que si vous redirigez la sortie standard dans un fichier directement dans la commande clingo pour plus de latitude, la sortie standard du travail sera alors vide ou presque vide...
