---
id: sppider_overview
title: SPPIDER Overview
---

SPPIDER is a web server that can predict sites of protein-protein interaction.

Unfortunately, it is only available as a web server. Forunately, that web server will email you results in plain text.

So, we can use `sppider_mech_JAS.pl` to automatically submit jobs to the web server, outlook to download emails as .csv, and then `parseSPPIDER2results.pl` to parse those exported emails into a reasonable format.
