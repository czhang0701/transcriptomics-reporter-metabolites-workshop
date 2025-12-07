# 3-Hour Workshop Schedule
## Transcriptomics to Reporter Metabolites

**Total Duration**: 3 hours (with breaks)
**Focus**: Get to reporter metabolite analysis
**Format**: Hands-on, code-along

---

## ⏰ Detailed Timeline

### Pre-Workshop (Students Complete at Home)
**Time**: 45-60 minutes (BEFORE workshop day)
- Install R and RStudio
- Install all packages
- Download materials
- **See**: `PRE_WORKSHOP_SETUP.md`

---

### Workshop Day Schedule

| Time | Duration | Activity | Module | Type |
|------|----------|----------|--------|------|
| 0:00-0:10 | 10 min | Welcome & Setup Check | - | Intro |
| 0:10-0:40 | 30 min | Differential Expression | 1 (condensed) | Code-along |
| 0:40-1:00 | 20 min | Gene Set Enrichment | 2 (condensed) | Code-along |
| **1:00-1:10** | **10 min** | **BREAK** | - | Break |
| 1:10-1:30 | 20 min | Co-expression Networks | 3 (condensed) | Code-along |
| 1:30-2:20 | 50 min | ⭐ Reporter Metabolites | 4 (full) | Code-along |
| 2:20-2:30 | 10 min | Interpretation & Discussion | - | Discussion |
| **2:30-2:40** | **10 min** | **BREAK** | - | Break |
| 2:40-3:00 | 20 min | Q&A, Next Steps, MOFA2 Preview | Bonus | Discussion |

**Total**: 3 hours (includes 20 min breaks)

---

## 📋 Module Details

### 0:00-0:10 | Welcome & Setup Check (10 min)

**Instructor Activities**:
- Welcome participants
- Overview of workshop goals
- Check everyone has R/RStudio working

**Student Activities**:
```r
# Everyone run this in RStudio:
getwd()  # Where am I?
setwd("path/to/workshop_materials")  # Adjust path!

# Verify packages
source("scripts/00_check_setup.R")
# Should see: ✓✓✓ ALL PACKAGES INSTALLED
```

**Learning Objectives**:
- Understand workshop goals
- Confirm setup working
- Know where materials are located

**Key Message**: *"Today we're going from 20,000 genes to identifying key metabolic reprogramming events"*

---

### 0:10-0:40 | Module 1: Differential Expression (30 min)

**Focus**: Essential DESeq2 workflow ONLY

**Topics Covered**:
1. Load data (5 min)
2. Create DESeq object (5 min)
3. Run DESeq2 (5 min)
4. Extract results (5 min)
5. Volcano plot (5 min)
6. Save results for next steps (5 min)

**Code Sections** (from Module1 tutorial):
- Steps 1-7 ONLY
- SKIP: QC plots, detailed exploration
- FOCUS: Get to differential expression results quickly

**Key Files Created**:
- `DESeq_output.txt` - for GSEA and reporter metabolites

**Learning Objectives**:
- Understand DESeq2 input requirements
- Run differential expression
- Interpret log2FoldChange and padj
- Identify significant genes

**Student Exercise** (5 min within the 30 min):
- How many significant genes (padj < 0.05)?
- Quick answer, move on

**Key Message**: *"We have 2,000 significant genes. What biological processes are affected?"*

---

### 0:40-1:00 | Module 2: Gene Set Enrichment (20 min)

**Focus**: Fast pathway enrichment to identify metabolic changes

**Topics Covered**:
1. Map Ensembl to gene symbols (3 min)
2. Load gene sets (2 min)
3. Run PIANO (5 min - while running, explain concepts)
4. View top pathways (5 min)
5. Identify metabolic pathways (5 min)

**Code Sections** (from Module2 tutorial):
- Steps 2-7 ONLY
- SKIP: Detailed visualization, consensus scoring
- FOCUS: Get enrichment results quickly

**Key Files Created**:
- `Piano_output.txt` - pathway enrichment results

**Learning Objectives**:
- Understand why GSEA is powerful
- Run PIANO enrichment
- Identify enriched metabolic pathways
- See which pathways up/down regulated

**Quick Discussion** (within the 20 min):
- "Glycolysis is up-regulated - what does this mean?"
- "TCA cycle changes - why important?"

**Key Message**: *"Metabolic pathways are reprogrammed. Which specific metabolites drive this?"*

---

### 1:00-1:10 | BREAK (10 min)

**Activities**:
- Stretch
- Coffee/water
- Check everyone's code is working
- Troubleshoot any errors

**Instructor**: Use this time to help anyone who fell behind

---

### 1:10-1:30 | Module 3: Co-expression Networks (20 min)

**Focus**: Quick network analysis to identify gene modules

**Topics Covered**:
1. Load FPKM data (3 min)
2. Calculate correlations (5 min)
3. FDR correction (3 min)
4. Filter significant (3 min)
5. Visualize network (3 min)
6. Export for Cytoscape (3 min)

**Code Sections** (from Module3 tutorial):
- Steps 2-7, 10 ONLY
- SKIP: Module detection, detailed network analysis
- FOCUS: Get correlation network for context

**Key Files Created**:
- `Coexpression_network.txt` - for Cytoscape (optional later)

**Learning Objectives**:
- Understand co-expression concept
- See which genes are correlated
- Identify metabolic gene clusters

**Quick Observation**:
- "Glycolysis genes are highly correlated"
- "This confirms coordinated regulation"

**Key Message**: *"Genes work together in metabolic networks. Now let's find the key metabolites."*

---

### 1:30-2:20 | ⭐ Module 4: Reporter Metabolites (50 min)

**Focus**: MAIN EVENT - Full reporter metabolite analysis

**THIS IS THE CORE OF THE WORKSHOP**

#### Part 1: Concept & Setup (10 min)

**Topics**:
- What are reporter metabolites? (5 min)
- How does the method work? (3 min)
- Why is this useful? (2 min)

**Slides/Discussion**: Instructor explains concept with examples

#### Part 2: Analysis (25 min)

**Topics**:
1. Load DESeq results (3 min)
2. Map to gene symbols (3 min)
3. Create metabolite-gene associations (5 min)
4. Run reporter analysis (8 min - explain while running)
5. Extract results (3 min)
6. Identify top reporters (3 min)

**Code Sections** (from Module4 tutorial):
- Steps 2-7 COMPLETELY
- Focus on PIANO reporter method
- Use example metabolite gene sets (faster than full model)

**Key Files Created**:
- `Reporter_Metabolites_output.txt` - THE MAIN RESULT!

#### Part 3: Interpretation (15 min)

**Activities**:
1. View top reporter metabolites (5 min)
2. Biological interpretation (5 min)
   - Lactate → Warburg effect
   - Glutamine → Glutamine addiction
   - ATP → Energy metabolism
3. Connect to earlier results (5 min)
   - Glycolysis enriched (Module 2) + Lactate reporter (Module 4)
   - TCA changes (Module 2) + TCA metabolites (Module 4)

**Group Discussion**:
- "What do these reporter metabolites tell us?"
- "How does this connect genes to metabolism?"
- "What would you measure experimentally?"

**Learning Objectives**:
- Understand reporter metabolite concept ✓
- Run reporter analysis ✓
- Interpret biological meaning ✓
- Connect to pathway enrichment ✓
- Identify metabolic reprogramming ✓

**Key Message**: *"We've identified key metabolic changes at the metabolite level - this is the mechanistic interpretation!"*

---

### 2:20-2:30 | Interpretation & Wrap-up (10 min)

**Activities**:
1. **Review what we learned** (3 min):
   - Genes → Pathways → Networks → Metabolites
   - Complete analysis pipeline
   - Biological interpretation

2. **Key Takeaways** (3 min):
   - Reporter metabolites link transcriptomics to metabolism
   - Provides mechanistic insights
   - Identifies therapeutic targets

3. **Real-world application** (4 min):
   - How to apply to your data
   - When to use reporter metabolites
   - Resources for genome-scale models

**Deliverables**: Students now have:
- ✓ Differential expression results
- ✓ Enriched pathways
- ✓ Co-expression network
- ✓ Reporter metabolites ⭐
- ✓ Biological interpretation

---

### 2:30-2:40 | BREAK (10 min)

**Activities**:
- Final stretch
- Prepare for Q&A
- Save all your work!

---

### 2:40-3:00 | Q&A & Next Steps (20 min)

#### Q&A Session (10 min)
- Answer questions from workshop
- Clarify concepts
- Troubleshoot issues

#### Next Steps (5 min)
**For continued learning**:
1. **Practice**: Use practice dataset
2. **Apply**: To your own data
3. **Explore**: Full tutorials (Modules 1-4 detailed versions)
4. **Advanced**: MOFA2 bonus modules (5-8)

**Resources provided**:
- Full tutorial PDFs
- Example scripts
- Solutions to exercises
- Contact information

#### MOFA2 Preview (5 min)
**Quick introduction to bonus content**:
- What is MOFA2?
- When to use it?
- How it extends reporter metabolite analysis
- Self-study materials available

**Key Message**: *"You now have the foundation. Continue learning with the full materials!"*

---

## 📊 What Students Will Learn

### By End of Workshop

**Conceptual Understanding**:
- ✓ Differential expression analysis workflow
- ✓ Gene set enrichment analysis (GSEA)
- ✓ Co-expression network concepts
- ✓ **Reporter metabolite analysis** ⭐
- ✓ Genome-scale metabolic models
- ✓ Transcriptomics to metabolism integration

**Practical Skills**:
- ✓ Run DESeq2 analysis
- ✓ Perform pathway enrichment with PIANO
- ✓ Build correlation networks
- ✓ **Identify reporter metabolites** ⭐
- ✓ Interpret metabolic reprogramming
- ✓ Generate publication-quality results

**Analytical Thinking**:
- ✓ Connect genes → pathways → metabolites → biology
- ✓ Interpret multi-level omics results
- ✓ Generate testable hypotheses
- ✓ Identify therapeutic targets

---

## 🎯 Workshop Success Criteria

**Minimum Success** (everyone should achieve):
- [ ] Completed Module 4 (Reporter Metabolites)
- [ ] Generated reporter metabolite rankings
- [ ] Understood biological interpretation
- [ ] Know how to apply to own data

**Ideal Success**:
- [ ] Completed all 4 modules
- [ ] Generated all output files
- [ ] Understood connections between modules
- [ ] Can explain results to others

**Bonus** (if time permits):
- [ ] Explored MOFA2 concepts
- [ ] Started practice dataset analysis
- [ ] Discussed own research applications

---

## 📚 Materials Organization

### Students Should Have Open

**During Workshop**:
1. **RStudio** - for coding
2. **Tutorial PDFs** - as reference (if desired)
3. **This schedule** - to track progress

**Files to Navigate To**:
```
workshop_materials/
├── data/              # All data files here
├── scripts/           # Helper scripts
├── tutorials/         # Full detailed tutorials (reference)
└── figures/           # Your plots will save here
```

---

## 👨‍🏫 Instructor Notes

### Pacing Strategy

**Strict Time Management**:
- Set timer for each module
- Move on even if some students lag (help during breaks)
- MUST reach Module 4 - this is the goal

**Code-Along Approach**:
- Type code live (don't copy-paste)
- Explain each line briefly
- Students type along
- Quick checks: "Does everyone see X?"

**Break Management**:
- Use breaks to help struggling students
- Keep breaks SHORT (10 min max)
- Resume promptly

### Adaptation Strategies

**If Running Ahead of Schedule**:
- Add more interpretation/discussion
- Show additional visualizations
- Answer more questions
- Preview MOFA2 in more detail

**If Running Behind**:
- Skip visualization steps
- Provide pre-made results
- Focus on concepts over code
- Extend Module 4 if needed

**Priority Order**:
1. Module 4 (reporter metabolites) - MUST COMPLETE
2. Module 2 (GSEA) - very important for context
3. Module 1 (DESeq2) - can use pre-made results
4. Module 3 (networks) - can skip if needed

---

## 💻 Technical Requirements

### During Workshop

**Each Student Needs**:
- Laptop with R/RStudio installed
- All packages installed (from PRE_WORKSHOP_SETUP.md)
- Workshop materials downloaded
- Stable internet (for troubleshooting)

**Instructor Needs**:
- Projector/screen sharing
- Backup materials on USB
- Pre-run code with outputs (for emergencies)
- List of common errors/solutions

---

## 🔧 Troubleshooting Plan

### Common Issues

**"Package X not found"**:
- Quick install during break
- Pair student with neighbor
- Use instructor's screen

**"Code doesn't work"**:
- Check working directory
- Check file paths
- Re-run from beginning

**"Results look different"**:
- Likely OK - biological data varies
- Check if error or just different numbers

**"Too fast!"**:
- Encourage pairing with faster student
- Provide pre-made results for that step
- Help during breaks

---

## 📧 Post-Workshop Follow-up

**Immediately After**:
- Share all materials via email/website
- Provide feedback survey
- Offer office hours for questions

**Within 1 Week**:
- Send solutions to exercises
- Share additional resources
- Invite to online discussion group

**Ongoing**:
- Email support for applying to own data
- Share updates to materials
- Connect students with each other

---

## ✅ Pre-Workshop Instructor Checklist

**1 Week Before**:
- [ ] Send PRE_WORKSHOP_SETUP.md to all students
- [ ] Test all code on Windows and Mac
- [ ] Prepare backup pre-run results
- [ ] Create Zoom/Teams meeting (if online)

**1 Day Before**:
- [ ] Email reminder with setup verification
- [ ] Prepare presentation slides
- [ ] Test screen sharing
- [ ] Print this schedule

**Workshop Day**:
- [ ] Arrive 15 minutes early
- [ ] Test projection/sharing
- [ ] Have USB backup of materials
- [ ] Start on time!

---

## 🎓 Learning Outcomes Assessment

### Formative (During Workshop)

**Quick Checks Every 20 Minutes**:
- "Raise hand if you see the plot"
- "Type 'yes' in chat if code ran"
- "What is the top enriched pathway?"

### Summative (End of Workshop)

**Final Questions** (during Q&A):
1. "What are reporter metabolites?"
2. "Why is lactate a reporter in cancer?"
3. "How does this connect to GSEA results?"

**Success Indicators**:
- Students can explain reporter metabolites
- Students generated all key output files
- Students understand biological interpretation

---

## 📖 Additional Resources for Students

**Provided in Workshop Materials**:
- Full detailed tutorials (Modules 1-4)
- Practice dataset and scripts
- Solutions to exercises
- Links to papers and databases

**Recommended Reading** (optional):
- Patil & Nielsen (2005) - Reporter metabolites paper
- DESeq2 vignette
- PIANO documentation

**For Continued Learning**:
- MOFA2 bonus modules (self-paced)
- Apply to own data
- Join Bioconductor community

---

## 🎯 Final Checklist for Students

**At End of Workshop, You Should Have**:

**Files Generated**:
- [ ] `DESeq_output.txt`
- [ ] `Piano_output.txt`
- [ ] `Coexpression_network.txt`
- [ ] `Reporter_Metabolites_output.txt` ⭐

**Knowledge Gained**:
- [ ] Can run differential expression
- [ ] Understand pathway enrichment
- [ ] Know what reporter metabolites are
- [ ] Can interpret metabolic reprogramming

**Next Steps**:
- [ ] Practice with practice dataset
- [ ] Work through full tutorials
- [ ] Apply to own data
- [ ] Explore MOFA2 (bonus)

---

**Workshop designed to guarantee completion of reporter metabolite analysis in 3 hours!** 🎯

**Questions?** See instructor or email: [contact email]
