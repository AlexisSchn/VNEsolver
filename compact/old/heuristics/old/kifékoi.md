A tester sur une instance coeur sur coeur et sur une instance subgraphs.


- heuristic_relax: 
basique, on prend les meilleurs x, on restreint dessus, et on balance à la résol.
Résultats navrants...

- heuristic_better:
On met les star csts. A verif qu'elles marchent toutes bieng...
Donne des résultats meilleurs.
Mais quand même naze sur des graphes peu structurés il semble...

- heuristique_better_progressive
on "branche" pour un noeud virtuel et on refait le PL...

- heuristic y:
pas encore fait. L'idée c'est d'enlever des y pareillement parce que ça peut être trop long à la résol.


Enfin bon de toute façon l'heuristique est bcp trop nulle de base pour que ce soit intéressant d'enlever des y en plus ptdr. Même en gardant 30 noeuds (sur 70 certes) je suis à l'ouest... Meme en diving... Enorme flop cette heuristique, qui montre que vine serait affreux pour ces instances.