# Documentation Updates Needed

## Module Reorganization Complete

**Main.R has been reorganized:**
- Module 1: DESeq2 (20 min)
- Module 2: PIANO/GSEA (15 min)
- **Module 3: Reporter Metabolites** ⭐ (60 min) - WAS Module 4
- **Module 4: Advanced Co-expression** (30 min) - REPLACED old basic Module 3

## Files Updated ✅
- [x] Main.R - Completely reorganized
- [x] README.md - Timeline and module table updated
- [x] scripts/coexpression_modules.R - New advanced functions

## Files Still Need Updates 📝

### High Priority:
- [ ] WORKSHOP_GUIDE_COMPLETE.md
  - Lines 270+: Module 1 (OK)
  - Lines 361+: Module 2 (OK)
  - Lines 436+: Module 3 currently says "Co-expression" → should be "Reporter Metabolites"
  - Need to add: Module 4 "Advanced Co-expression"
  - Need to remove: Old basic co-expression content

### Tutorial Files:
- [ ] tutorials/Module3_Coexpression.md → Rename to Module4_Coexpression.md
- [ ] tutorials/Module4_Reporter Metabolites.md → Rename to Module3_ReporterMetabolites.md
- [ ] Update internal references in tutorials

### Medium Priority:
- [ ] CLAUDE.md (in parent directory) - Update module descriptions
- [ ] deprecated/ files - Update if needed

### Low Priority:
- [ ] Solution files - Update module numbers if they exist

## Quick Fix Commands

```bash
# Rename tutorial files
mv tutorials/Module3_Coexpression.md tutorials/Module4_Coexpression.md
mv tutorials/Module4_ReporterMetabolites.md tutorials/Module3_ReporterMetabolites.md

# Update internal references
sed -i 's/Module 3: Co-expression/Module 4: Advanced Co-expression/g' WORKSHOP_GUIDE_COMPLETE.md
sed -i 's/Module 4: Reporter/Module 3: Reporter/g' WORKSHOP_GUIDE_COMPLETE.md
```

## Notes
- Main functionality is correct (Main.R works)
- Documentation just needs to catch up with code
- Not breaking - just inconsistent numbering in docs
