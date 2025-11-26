# Leptospirosis outbreak thresholds

**General disclaimer** This repository was created for use by CDC programs to collaborate on public health related projects in support of the [CDC mission](https://www.cdc.gov/about/cdc/#cdc_about_cio_mission-our-mission).  GitHub is not hosted by the CDC, but is a third party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise. 


## Overview

This repository provides code and supporting documentation for deriving weekly leptospirosis outbreak thresholds in Puerto Rico using historical case data and a statistical modeling framework, with an adjustment to allow for the change in reporting definition in 2022. The approach informs real-time public health decision-making in collaboration with the Puerto Rico Department of Health. The approach is a modified version of a previously cited method (https://bmjopen.bmj.com/content/15/9/e106182) that allows for estimation with lower case counts.

### Background

Leptospirosis remains a recurring public health concern in Puerto Rico and other endemic regions, particularly following hurricanes and flooding events when environmental conditions favor increased exposure. Here, I have adapted disease outbreak threshold methods previously published in BMJ Open (https://bmjopen.bmj.com/content/15/9/e106182) to allow for low weekly case counts. This method uses a weekly **intercept-only negative binomial regression model**, fit to historical data using time-series bootstrapping to estimate threshold values. These thresholds enable objective, data-driven identification of outbreak weeks and timely alerts that can be incoporated into reoutine surveillance reports and used to monitor leptospirosis activity.

### Repository Contents


* **`/code/threshold_model_lepto.R`** - (.R file) Primary R script for running thresholds and plotting figures. The bullets below are done at both the island-wide and health region level
  * Reads in the data, cleans and formats it on a weekly scale.
  * Draws time-series bootstrap estimates to supplement the dataset.
  * Fits the intercept-only negative binomial model using the bootstrapped data.
  * Calculates percentile-based thresholds (75th, 80th, 85th, 90th, and 95th).
  * Applies a range of consecutive week rules (2+, 3+, and 4+).
  * Generates graphs.
* **`/.github/`** - (.md files) issue and pull request templates

### Usage Instructions

1. **Locate data**: Place the case data files in the `/data` folder, as an R-readable file (e.g., CSV, RDS). Alternatively, note the file path where they are saved to be read in later.
2. **Run the model**: Execute `threshold_model_lepto.R`. This will compute thresholds and generate key figures.


These examples demonstrate practical impact and assist in bridging surveillance data with actionable public health responses.

### Citation

If you use these materials or replicate the modeling approach, please cite this code and:
Thayer MB, Marzan-Rodriguez M, Torres Aponte J, et al. Dengue epidemic alert thresholds for surveillance and decision-making in Puerto Rico: development and prospective application of an early warning system using routine surveillance data. BMJ Open 2025; 15:e106182. doi: 10.1136/bmjopen-2025-106182.



  
## Public Domain Standard Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.


The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](DISCLAIMER.md)
and [Code of Conduct](code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template) for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/main/CONTRIBUTING.md), [public domain notices and disclaimers](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md), and [code of conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md).
