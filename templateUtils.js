// Template Utils

var templateUtils = templateUtils || {};

templateUtils.defaults = _.partialRight(_.assign, function (value, other) {
  return _.isUndefined(value) ? other : value;
});

templateUtils.labelData = function (input, layoutSection) {
  var obj = {};
  obj.id = layoutSection.id;
  obj.label = layoutSection.label;
  obj.data = input[layoutSection.id];
  return obj;
};

templateUtils.combineLayoutWithData = function (input, layoutSection) {
  var obj = {};
  var key = layoutSection.id;
  var data = input[key];
  var curried = _.curry(templateUtils.labelData)(data);
  var labeledData = _.map(layoutSection.items, curried);
  obj[key] = labeledData;
  return obj;
};

templateUtils.labelInput = function (layout, input) {
  var curried = _.curry(templateUtils.combineLayoutWithData)(input);
  // label all inputs except those with variant context
  // var sections = _.map(_.filter(layout, function (section) { return section.context !== 'variant'; }), curried);
  var sections = _.map(_.filter(layout, function (section) { return ['variant', 'cnv'].indexOf(section.context) < 0; }), curried);

  var labeledSections = _.reduce(sections, _.extend, {});
  // now add back the variant arrays
  delete labeledSections['cnvs'];
  labeledSections = templateUtils.defaults(labeledSections, input);
  return labeledSections;
};

templateUtils.chunkObject = function (input, size) {
  // pretty certain this is a bug in lodash, 
  // but we need at separate scope for the _.object method... i'll dig into it and submit a PR
  var closureToObject = function (arr) {
    return _.object(arr);
  };
  return _(input).pairs().chunk(size).map(closureToObject).value();
};

templateUtils.divideObject = function (obj, prop, pieces) {
  var subObj = obj[prop];
  var size = Math.ceil(Object.keys(subObj).length / pieces);
  var chunks = templateUtils.chunkObject(subObj, size);
  // delete obj[prop];
  _.forEach(chunks, function (chunk, idx) {
    obj[prop + '_' + idx.toString()] = chunk;
  });
  return obj;
};

templateUtils.formatFrequency = function (freq) {
  freq *= 100;
  if (freq < 1) {
    freq = _.round(freq, 4);
  } else if (freq < 10) {
    freq = _.round(freq, 2);
  } else if (freq < 90) {
    freq = _.round(freq);
  } else if (freq < 99) {
    freq = _.round(freq, 2);
  } else {
    freq = _.round(freq, 4);
  }
  return freq + "%";
};

templateUtils.arrayToString = function (arr) {
  var size = arr.length;
  if (size === 0) {
    return '';
  } else if (size === 1) {
    return arr[0];
  } else if (size === 2) {
    return arr.join(' and ');
  } else {
    arr[size - 1] = 'and ' + arr[size - 1];
    return arr.join(', ');
  }
};

// templateUtils.buildGeneDescriptionHash = function(variants) {
//   // predicate
//   var omimGenesArrPicker = function(value, key){ 
//       return _.startsWith(key, "OMIM-Genes") && Array.isArray(value); };

//   var desc = {};
//   _.forEach(variants, function(variant){
//         primaryGeneKeys = _.chain(variant).pick(omimGenesArrPicker).values().flatten()
//           .pluck('GeneName').value()[0].split(',');
//         altGeneKeys = _.chain(variant).pick(omimGenesArrPicker).values().flatten()
//             .pluck('AlternativeGeneSymbols').value()[0].split(',');
//         geneKeys = _.union(primaryGeneKeys, altGeneKeys);
//         geneDesc = _.chain(variant).pick(omimGenesArrPicker).values().flatten().value()[0];
//         comboHash = _.zipObject(geneKeys, _.fill(Array(geneKeys.length), geneDesc));
//         desc = _.extend(desc, comboHash);
//   });
//   return desc;
// };


templateUtils.buildPhenotypeHash = function (variants) {
  // predicate
  var omimPhenotypeArrPicker = function (value, key) {
    return _.startsWith(key, "OMIM-Phenotypes") && Array.isArray(value);
  };

  var phash = {};
  _.forEach(variants, function (variant) {
    pheno = _.chain(variant).pick(omimPhenotypeArrPicker).values().flatten().groupBy('GeneNames').value();
    phash = _.extend(phash, pheno);
  });
  return phash;
};

templateUtils.buildVariantHash = function (variants) {
  //predicate
  var omimVarArrPicker = function (value, key) {
    return _.startsWith(key, "OMIM-Variants") && Array.isArray(value);
  };

  var varHash = {};
  _.forEach(variants, function (variant) {
    varDesc = _.chain(variant).pick(omimVarArrPicker).values().flatten().indexBy('GeneName').value();
    varDesc = _.extend(varHash, varDesc);
  });
  return varHash;
};

templateUtils.buildOncoMDDrugArray = function (variant) {
  //predicate
  var oncoMDDrugArrPicker = function (value, key) {
    return _.startsWith(key, "VariantDrugs-OncoMD") && Array.isArray(value);
  };

  var drugArr = _.chain(variant).pick(oncoMDDrugArrPicker).values().flatten().value();
  return drugArr;
};

templateUtils.buildOncoMDTrialsArray = function (variant) {
  //predicate
  var oncoMDTrialsArrPicker = function (value, key) {
    return _.startsWith(key, "ClinicalTrials-OncoMD") && Array.isArray(value);
  };

  var trialsArr = _.chain(variant).pick(oncoMDTrialsArrPicker).values().flatten().value();
  return trialsArr;
};

templateUtils.omimString = function (variant) {
  var omimStr = '';

  // var geneDescHash = templateUtils.buildGeneDescriptionHash([variant]);
  var phash = templateUtils.buildPhenotypeHash([variant]);
  var vhash = templateUtils.buildVariantHash([variant]);

  _.forEach(variant.transcript.GeneNames, function (geneName) {
    omimStr += _.chain(vhash).get([geneName, 'Description'], '').value();

    // var geneDesc = _.get(geneDescHash, geneName, {});
    // omimStr += _.get(geneDesc, 'Description', '');

    // var disorderInteritanceArr = _.chain(geneDesc).get('DisordersInheritance').unique().value();
    // if(disorderInteritanceArr.length > 0){
    //   omimStr += '<p>';
    //   omimStr += 'This gene has been observed to exhibit ';
    //   omimStr += templateUtils.arrayToString(disorderInteritanceArr);
    //   omimStr += disorderInteritanceArr > 1 ? ' inheritance patterns.' : ' inheritance pattern.';
    //   omimStr += '</p>';
    // }

    // var disorders = _.get(geneDesc, 'Disorders', []);
    // if(disorders.length > 0){
    //   omimStr += '<p>';
    //   omimStr += 'It has been associated with ' + 
    //     templateUtils.arrayToString(disorders) + '.';
    //   omimStr += '</p>';
    // }
    var phenoDescArr = _.chain(phash).get(geneName).pluck('Description').filter(null).value();
    omimStr += phenoDescArr.join('\n');
  }); // foreach gene
  return omimStr;
};

templateUtils.omimGeneString = function (variant) {
  var omimStr = '';

  var geneDescHash = templateUtils.buildGeneDescriptionHash([variant]);
  var phash = templateUtils.buildPhenotypeHash([variant]);

  _.forEach(variant.transcript.GeneNames, function (geneName) {
    var geneDesc = _.get(geneDescHash, geneName, {});
    omimStr += _.get(geneDesc, 'Description', '');

    var disorderInteritanceArr = _.chain(geneDesc).get('DisordersInheritance').unique().value();
    if (disorderInteritanceArr.length > 0) {
      omimStr += '<p>';
      omimStr += 'This gene has been observed to exhibit ';
      omimStr += templateUtils.arrayToString(disorderInteritanceArr);
      omimStr += disorderInteritanceArr > 1 ? ' inheritance patterns.' : ' inheritance pattern.';
      omimStr += '</p>';
    }

    var disorders = _.get(geneDesc, 'Disorders', []);
    if (disorders.length > 0) {
      omimStr += '<p>';
      omimStr += 'It has been associated with ' +
        templateUtils.arrayToString(disorders) + '.';
      omimStr += '</p>';
    }
    var phenoDescArr = _.chain(phash).get(geneName).pluck('Description').filter(null).value();
    omimStr += phenoDescArr.join('\n');
  }); // foreach gene

  return omimStr;
};

templateUtils.soString = function (variant) {
  var soStr = '';
  var soStrArr = [];
  soStrArr = _.chain(variant).get('transcript.SequenceOntologyClinicallyRelevant', []).map(function (so) {
    return _.chain(so).split('_').map(_.capitalize).value().join(' ');
  }).value();
  var tx = _.get(variant, 'transcript.TranscriptNameClinicallyRelevant', []);
  var genes = _.get(variant, 'transcript.GeneNames', []);
  if ((genes) && (genes.length === 0)) {
    soStr = 'This is an Intergenic Variant.';
  } else if (soStrArr.length === 1) {
    soStr = 'This is a ' + soStrArr[0] + ' located in the ' + genes[0] + ' gene.';
  } else if (soStrArr.length > 1) {
    soStr = 'This is a ';
    _.forEach(soStrArr, function (so, idx) {
      soStr += so + ' (' + tx[idx] + ')';
      if (idx == soStrArr.length - 1)
        soStr += '. ';
      else if (idx == soStrArr.length - 2)
        soStr += ' and ';
      else
        soStr += ', ';
    });
    soStr += "It is located in the ";
    soStr += '<i>' + templateUtils.arrayToString(genes) + '</i>';
    soStr += genes.length > 1 ? ' genes.' : 'gene.';
  }
  soStr += '\n';
  return soStr;
};

templateUtils.soVariantString = function (variant) {
  var soStr = '';
  var soStrArr = _.chain(variant).get('transcript.SequenceOntologyClinicallyRelevant', []).map(function (so) {
    return _.chain(so).split('_').map(_.capitalize).value().join(' ');
  }).value();
  var tx = _.get(variant, 'transcript.TranscriptNameClinicallyRelevant', []);
  var genes = _.get(variant, 'transcript.GeneNames', []);
  if ((genes) && (genes.length === 0)) {
    soStr = 'This is an Intergenic Variant.';
  } else if (soStrArr.length === 1) {
    soStr = 'This is a ' + soStrArr[0] + '.';
  } else if (soStrArr.length > 1) {
    soStr = 'This is a ';
    _.forEach(soStrArr, function (so, idx) {
      soStr += so + ' (' + tx[idx] + ')';
      if (idx == soStrArr.length - 1)
        soStr += '. ';
      else if (idx == soStrArr.length - 2)
        soStr += ' and ';
      else
        soStr += ', ';
    });
  }
  soStr += '\n';
  return soStr;
};

templateUtils.varFreqNHLBIString = function (variant) {
  var freqStr = '';
  var freq = _.get(variant, 'NHLBI.AllAAF');
  if (freq) {
    freqStr = "It appears with a frequency of " +
      templateUtils.formatFrequency(freq) +
      " in the NHLIBI 6500 Exomes Project.";
  }
  return freqStr;
};

templateUtils.varFreq1KGString = function (variant) {
  var freqStr = '';
  var freq = _.get(variant, '1kG.AllIndivFreq');
  if (freq) {
    freqStr = "It appears with a frequency of " +
      templateUtils.formatFrequency(freq) +
      " in the 1000 Genomes Project.";
  }
  return freqStr;
};

templateUtils.varFreqExacString = function (variant) {
  var freqStr = '';
  var freqArr = _.get(variant, 'ExAC-Variant Frequencies.AdjustedAltAlleleFreqAF_Adj', []);
  var freqs = templateUtils.arrayToString(_.map(freqArr, templateUtils.formatFrequency));
  if (freqArr.length == 1) {
    freqStr = "It appears with a frequency of " + freqs +
      " in the ExAC Variant Catalog.";
  } else if (freqArr.length > 1) {
    freqStr = "It appears with frequencies of " + freqs +
      " in the ExAC Variant Catalog.";
  }
  return freqStr;
};

templateUtils.variantTitle = function (variant, input) {
  var titleTemplate;
  var titles = [];
  var genes = _.get(variant, 'transcript.GeneNames', []);
  _.forEach(genes, function (gene, idx) {
    var txPair = {
      'chr': variant.site.Chr,
      'start': variant.site.Start,
      'refalt': variant.site.RefAlt,
      'gene': gene,
      'exon': _.get(variant, ['transcript', 'ExonNumberClinicallyRelevant', idx], ''),
      'hgvsC': _.get(variant, ['transcript', 'HGVScClinicallyRelevant', idx], ''),
      'hgvsP': _.get(variant, ['transcript', 'HGVSpClinicallyRelevant', idx], ''),
      'pathogenicity': input.classification,
      'sequenceOntology': _.chain(variant).get(['transcript', 'SequenceOntologyClinicallyRelevant', idx], '').split('_').map(_.capitalize).value().join(' ')
    };
    if (txPair.hgvsP !== '?') {
      titleTemplate = _.template('{{hgvsP}} in Exon {{exon}} of <i>{{gene}}</i> ({{hgvsC}}) {{pathogenicity}}');
    } else if (txPair.gene !== '?') {
      titleTemplate = _.template('{{sequenceOntology}} in <i>{{gene}}</i> ({{hgvsC}}) {{pathogenicity}}');
    } else {
      titleTemplate = _.template('{{sequenceOntology}} (Chr{{chr}}:{{start}} {{refalt}}) {{pathogenicity}}');
    }
    titles.push('<h4>' + titleTemplate(txPair) + '</h4>');
  });
  if ((genes) && (genes.length === 0)) {
    var hgvsG = _.get(variant, ['transcript_2', '0', 'HGVSg'], null);
    var rawNotation = variant.site.Chr + ':' + (variant.site.Start + 1) + ' ' + variant.site.RefAlt;
    var txPair = {
      'chr': variant.site.Chr,
      'start': variant.site.Start,
      'refalt': variant.site.RefAlt,
      'hgvs': hgvsG ? hgvsG : rawNotation,
      'pathogenicity': input.classification
    };
    titleTemplate = _.template('Intergenic Variant ({{hgvs}}) {{pathogenicity}}');
    titles.push('<h4>' + titleTemplate(txPair) + '</h4>');
  }
  return titles.join('\n');
};

templateUtils.getReferences = function (records, variants) {
  // Helper Functions
  var omimArrPicker = function (value, key) {
    return (_.startsWith(key, "OMIM-Genes") || _.startsWith(key, "OMIM-Phenotypes") ||
      _.startsWith(key, "OMIM-Variants")) && Array.isArray(value);
  };
  var collapseOnOmim = function (arr) {
    return _.values(_.pick(arr, omimArrPicker));
  };

  // regexs used in url parsing
  var urlRe = /href="([^\'\"]+)/g;
  var stripNumber = _.partial(String.prototype.replace, /^\d+\.\s/, '');

  // First, build up array of urls that are used in the interpretation sections
  var usedUrls = [];
  var interpConcat = _.chain(records).pluck('interpretation').value()
    .concat(_.chain(records).pluck('geneInterpretation').value())
    .concat(_.chain(records).pluck('variantInterpretation').value()).join();
  while (match = urlRe.exec(interpConcat)) {
    url = match[1].replace('href="', 'target="_blank" href="');
    usedUrls.push(url);
  }

  usedUrls = _.unique(usedUrls);
  // Now, make a hash of urls to refereces
  var referencesArr = _.chain(variants).values().map(collapseOnOmim)
    .flattenDeep().pluck('References').flattenDeep().omit(_.isNull).values().value();

  var refHash = {};
  _.forEach(referencesArr, function (ref) {
    while (match = urlRe.exec(ref)) {
      ref = ref.replace('<a href="', '<a target="_blank" href="');
      refHash[match[1]] = ref;
    }
  });

  // pick values in hash that are in the usedUrls arr, remove nulls, uniquify, strip id, and sort
  var usedReferences = _.chain(refHash).pick(usedUrls).values()
    .reject(_.isNull).reject(_.isUndefined).invoke(stripNumber).unique().sort().value();
  return usedReferences;
};

// CNV Template Utils

templateUtils.cnvTitle = function (region, input) {
  var titleTemplate = _.template('{{span}}bp {{type}} covering {{targets}} targets (Chr{{chr}}:{{start}}-{{stop}}) ');
  return titleTemplate({
    'chr': region.cnvcaller.Chr,
    'start': region.cnvcaller.Start,
    'stop': region.cnvcaller.Stop,
    'type': region.cnvcaller.Type,
    'targets': region.cnvcaller.Targets,
    'span': region.cnvcaller.Span
  });
};

templateUtils.buildGeneRegionDescriptionHash = function (region) {
  // predicate
  var fields = ['GeneNames', 'CDSCoveredClinicallyRelevant', 'CoveredClinicallyRelevant',
    'OverlappingExonsClinicallyRelevant', 'TranscriptNameClinicallyRelevant'];

  var genesArrPicker = function (value, key) {
    if (!_.startsWith(key, "overlapgene"))
      return false;
    return _.all(_.map(fields, _.curry(_.has)(value)));
  };

  var geneDesk = function (gene) {
    return _.zipObject(gene.GeneName, _.fill(Array(geneKeys.length), gene));
  };


  var genes = _.map(_.chain(region).pick(genesArrPicker).value(), function (overlapGene) {
    var geneList = [];
    for (i = 0; i < overlapGene.Genes; i++) {
      geneList.push(
        {
          'gene': overlapGene.GeneNames[i],
          'cds': overlapGene.CDSCoveredClinicallyRelevant[i],
          'covered': overlapGene.CoveredClinicallyRelevant[i],
          'overlapExons': overlapGene.OverlappingExonsClinicallyRelevant[i],
          'tx': overlapGene.TranscriptNameClinicallyRelevant[i]
        }
      );
    }
    return geneList;
  });

  return _.flatten(genes);
};


templateUtils.buildOmimGeneRegionDescriptionHash = function (regions) {
  // predicate
  var omimGenesArrPicker = function (value, key) {
    return _.startsWith(key, "overlap_OMIM-Genes") && Array.isArray(value);
  };

  var geneHash = function (gene) {
    var geneKeys = _.chain(gene).get('AlternativeGeneSymbols').split(',').union([gene.GeneName]).value();
    return _.zipObject(geneKeys, _.fill(Array(geneKeys.length), gene));
  };

  var desc = {};
  _.forEach(regions, function (region) {
    _.forEach(_.chain(region).pick(omimGenesArrPicker).values().flatten().map(geneHash).value(), function (geneDesc) {
      desc = _.extend(desc, geneDesc);
    });
  });

  return desc;
};

templateUtils.buildRegionPhenotypeHash = function (regions) {
  // predicate
  var omimPhenotypeArrPicker = function (value, key) {
    return _.startsWith(key, "overlap_OMIM-Phenotypes") && Array.isArray(value);
  };

  var phash = {};
  _.forEach(regions, function (region) {
    pheno = _.chain(region).pick(omimPhenotypeArrPicker).values().flatten().groupBy('GeneName').value();
    phash = _.extend(phash, pheno);
  });
  return phash;
};

templateUtils.buildRegionHash = function (regions) {
  //predicate
  var omimRegionArrPicker = function (value, key) {
    return key.indexOf("OMIM-Regions") !== -1 && Array.isArray(value);
  };


  // todo - decide how to organize these
  var regionHash = {};
  _.forEach(regions, function (region) {
    varDesc = _.chain(region).pick(omimRegionArrPicker).values().flatten().indexBy('GeneName').value();
    varDesc = _.extend(varHash, varDesc);
  });
  return varHash;
};

templateUtils.omimGeneRegionSummaries = function (region) {
  var omimGenesArrPicker = function (value, key) {
    return _.startsWith(key, "overlap_OMIM-Genes") && Array.isArray(value);
  };

  var geneHash = function (gene) {
    var geneKeys = _.chain(gene).get('AlternativeGeneSymbols').split(',').union([gene.GeneName]).value();

    var disorderCount = 0;
    if (gene.Disorders)
      disorderCount = gene.Disorders.length;

    var inheritanceType = "";
    if (gene.DisordersInheritance)
      inheritanceType = gene.DisordersInheritance;

    var summary = "Associated with " +
      disorderCount +
      (disorderCount > 0 ? " disorders" : " disorders") +
      " in OMIM with inheritance type \"" +
      inheritanceType + "\"";

    var geneList = [];
    _.forEach(geneKeys, function (gene) {
      geneList.push({ 'gene': gene, 'summary': summary });
    });

    return geneList;
  };

  return _.chain(region).pick(omimGenesArrPicker).values().flatten().map(geneHash).flatten().value();
};

templateUtils.omimGeneRegionPhenotypeStrings = function (region) {
  var geneDescHash = templateUtils.buildOmimGeneRegionDescriptionHash([region]);
  var phash = templateUtils.buildRegionPhenotypeHash([region]);

  var omimStrs = [];
  _.forEach(region.overlapgene.GeneNames, function (geneName) {

    var gene = _.get(geneDescHash, geneName, {});

    var geneDesc = _.get(gene, "Description");

    var title = _.chain(phash).get(geneName).pluck('Title').flatten().value();

    var omimStr = "<p>" + title + "</p>";

    if (_.has(geneDesc, 'Description'))
      omimStr += "<p>" + _.get(geneDesc, 'Description', '') + "</p>";

    var disorderInteritanceArr = _.chain(phash).get(geneName).pluck('Inheritance').flatten().unique().value();
    if (disorderInteritanceArr.length > 0) {
      omimStr += '<p>';
      omimStr += 'This gene phenotype has been observed to exhibit ';
      omimStr += templateUtils.arrayToString(disorderInteritanceArr);
      omimStr += disorderInteritanceArr.length > 1 ? ' inheritance patterns.' : ' inheritance pattern.';
      omimStr += '</p>';
    }

    var disorders = _.get(geneDesc, 'Disorders', []);
    if (disorders.length > 0) {
      omimStr += '<p>';
      omimStr += 'It has been associated with ' +
        templateUtils.arrayToString(disorders) + '.';
      omimStr += '</p>';
    }

    var phenoDescArr = _.chain(phash).get(geneName).pluck('Description').filter(null).unique().value();
    omimStr += phenoDescArr.join('\n');

    omimStrs.push({ 'gene': geneName, 'omimPhenotype': omimStr });
  });
  return omimStrs;
};

templateUtils.dosageClinGenString = function (cnv) {
  var cnvType = _.chain(cnv).get("Current").get("CNVState").value();
  if (cnvType == null || cnvType == "CN Loh")
    return "";

  var sensitivityString = "HaploinsufficiencyDescription";
  if (cnvType == "Duplicate")
    sensitivityString = "TriplosensitivityDescription";

  var clinGenArrPicker = function (value, key) {
    return key.indexOf("ClinGen") !== -1 && Array.isArray(value);
  };

  var filterType = function (clinGenRecord) {
    return ["No evidence available", "Not yet evaluated"].indexOf(
      _.get(sensitivityString, clinGenRecord, "Not yet evaluated")) != -1;
  };

  var dosageTemplate = _.template("<p>Observerd to exhibit '{{sensitivity}}' dosage sensitivity</p>");

  var dosageDescArr = _.chain(cnv).pick(clinGenArrPicker).values().flatten().filter(filterType).map(
    function (record) {
      return {
        'gene': _.get(record, "GeneSymbol"),
        'dosage': dosageTemplate({ 'sensitivity': _.get(record, sensitivityString) })
      };
    }).unique().value();

  return dosageDescArr;
};


