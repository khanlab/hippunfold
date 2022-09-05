# Pipeline Details

The HippUnfold workflow is generated based on the input data (e.g. whether there are multiple T2w images or a single T2w image), what modality is used (e.g. `modality=T1w` or `modality=T2w`), and what optional arguments are specified (e.g. `--t1-reg-template`). Below is a simplified rule graph of the **T1w** workflow (click on the image to enlarge):

<img src="../../hippunfold/dags/svg/T1w_rulegraph.svg"  width="1500px">

Each rounded rectangle in this diagram represents a *rule*, that is, some code or script that produces an output, and the arrows represent file inputs and outputs to these rules. Note that the `all` rule is special in that it is the target rule for the workflow (i.e. all the final output files of the workflow are inputs to this rule. The workflow diagram is also organized into groups of rules, which are defined by the names of the rule files, which can be found in the [rules sub-folder](http://github.com/khanlab/hippunfold/tree/master/hippunfold/workflow/rules)  in the workflow source. For example, the [preproc_t1](http://github.com/khanlab/hippunfold/tree/master/hippunfold/workflow/rules/preproc_t1.smk) file contains the rules related to pre-processing the T1w images, and these are grouped together in the above diagram by a blue rectangle labelled `preproc_t1`. 



