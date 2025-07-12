-- phpMyAdmin SQL Dump
-- version 5.2.1
-- https://www.phpmyadmin.net/
--
-- Host: 127.0.0.1
-- Generation Time: Jul 12, 2025 at 01:58 AM
-- Server version: 10.4.32-MariaDB
-- PHP Version: 8.2.12

SET SQL_MODE = "NO_AUTO_VALUE_ON_ZERO";
START TRANSACTION;
SET time_zone = "+00:00";


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8mb4 */;

--
-- Database: `subchems_search`
--

-- --------------------------------------------------------

--
-- Table structure for table `analysis_cache`
--

CREATE TABLE `analysis_cache` (
  `id` int(11) NOT NULL,
  `chemical_id` int(11) NOT NULL,
  `analysis_type` varchar(100) NOT NULL,
  `raw_data` longtext CHARACTER SET utf8mb4 COLLATE utf8mb4_bin DEFAULT NULL CHECK (json_valid(`raw_data`)),
  `analysis_date` timestamp NOT NULL DEFAULT current_timestamp(),
  `source` varchar(100) DEFAULT NULL,
  `version` varchar(50) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

--
-- Dumping data for table `analysis_cache`
--

INSERT INTO `analysis_cache` (`id`, `chemical_id`, `analysis_type`, `raw_data`, `analysis_date`, `source`, `version`) VALUES
(1, 6, 'comprehensive_analysis', '{\"basic_properties\":{\"smiles_canonical\":\"CCO\",\"molecular_formula\":\"C2H6O\",\"inchi\":\"InChI=1S\\/C2H6O\\/c1-2-3\\/h3H,2H2,1H3\",\"inchi_key\":\"LFQSCWFLJHTTHZ-UHFFFAOYSA-N\"},\"molecular_properties\":{\"molecular_weight\":46.069,\"exact_mass\":46.04186481,\"monoisotopic_mass\":46.04186481,\"xlogp3\":-0.0014,\"clogp\":-0.0014,\"hbd_count\":1,\"hba_count\":1,\"rotatable_bond_count\":0,\"heavy_atom_count\":3,\"aromatic_ring_count\":0,\"aliphatic_ring_count\":0,\"tpsa\":20.23,\"formal_charge\":0,\"bertz_complexity\":2.755,\"defined_atom_stereocenter_count\":0,\"undefined_atom_stereocenter_count\":0,\"covalently_bonded_unit_count\":1,\"compound_is_canonicalized\":true,\"has_stereochemistry\":false,\"radical_electron_count\":0,\"isotope_atom_count\":0,\"total_charge\":0},\"drug_properties\":{\"lipinski_violations\":0,\"lipinski_mw_violation\":false,\"lipinski_logp_violation\":false,\"lipinski_hbd_violation\":false,\"lipinski_hba_violation\":false,\"veber_violations\":0,\"veber_rotatable_bonds_violation\":false,\"veber_tpsa_violation\":false,\"qed_score\":0.4068,\"drug_like\":true,\"oral_bioavailability_likely\":true},\"physical_properties\":{\"estimated_logS\":0.04,\"estimated_melting_point\":null,\"estimated_boiling_point\":null},\"additional_descriptors\":{\"kappa1\":2.96,\"kappa2\":1.96,\"kappa3\":1.96,\"chi0v\":2.1543,\"chi1v\":1.0233,\"chi2v\":0.3162,\"ring_count\":0,\"saturated_ring_count\":0,\"heterocycle_count\":0,\"carbon_count\":2,\"nitrogen_count\":0,\"oxygen_count\":1,\"sulfur_count\":0,\"phosphorus_count\":0,\"halogen_count\":0},\"analysis_metadata\":{\"rdkit_version\":\"2025.03.3\",\"canonical_smiles\":\"CCO\",\"analysis_success\":true},\"success\":true}', '2025-07-11 20:24:12', 'rdkit', '2025.03.3'),
(2, 7, 'comprehensive_analysis', '{\"basic_properties\":{\"smiles_canonical\":\"O\",\"molecular_formula\":\"H2O\",\"inchi\":\"InChI=1S\\/H2O\\/h1H2\",\"inchi_key\":\"XLYOFNOQVPJJNP-UHFFFAOYSA-N\"},\"molecular_properties\":{\"molecular_weight\":18.015,\"exact_mass\":18.01056468,\"monoisotopic_mass\":18.01056468,\"xlogp3\":-0.8247,\"clogp\":-0.8247,\"hbd_count\":0,\"hba_count\":0,\"rotatable_bond_count\":0,\"heavy_atom_count\":1,\"aromatic_ring_count\":0,\"aliphatic_ring_count\":0,\"tpsa\":31.5,\"formal_charge\":0,\"bertz_complexity\":0,\"defined_atom_stereocenter_count\":0,\"undefined_atom_stereocenter_count\":0,\"covalently_bonded_unit_count\":1,\"compound_is_canonicalized\":true,\"has_stereochemistry\":false,\"radical_electron_count\":0,\"isotope_atom_count\":0,\"total_charge\":0},\"drug_properties\":{\"lipinski_violations\":0,\"lipinski_mw_violation\":false,\"lipinski_logp_violation\":false,\"lipinski_hbd_violation\":false,\"lipinski_hba_violation\":false,\"veber_violations\":0,\"veber_rotatable_bonds_violation\":false,\"veber_tpsa_violation\":false,\"qed_score\":0.3277,\"drug_like\":true,\"oral_bioavailability_likely\":true},\"physical_properties\":{\"estimated_logS\":0.732,\"estimated_melting_point\":null,\"estimated_boiling_point\":null},\"additional_descriptors\":{\"kappa1\":0.96,\"kappa2\":-27.04,\"kappa3\":-104.04,\"chi0v\":0.5,\"chi1v\":0,\"chi2v\":0,\"ring_count\":0,\"saturated_ring_count\":0,\"heterocycle_count\":0,\"carbon_count\":0,\"nitrogen_count\":0,\"oxygen_count\":1,\"sulfur_count\":0,\"phosphorus_count\":0,\"halogen_count\":0},\"analysis_metadata\":{\"rdkit_version\":\"2025.03.3\",\"canonical_smiles\":\"O\",\"analysis_success\":true},\"success\":true}', '2025-07-11 20:25:19', 'rdkit', '2025.03.3'),
(3, 8, 'comprehensive_analysis', '{\"basic_properties\":{\"smiles_canonical\":\"C=CCN(CC=C)CCc1c[nH]c2ccc(OC)cc12\",\"molecular_formula\":\"C17H22N2O\",\"inchi\":\"InChI=1S\\/C17H22N2O\\/c1-4-9-19(10-5-2)11-8-14-13-18-17-7-6-15(20-3)12-16(14)17\\/h4-7,12-13,18H,1-2,8-11H2,3H3\",\"inchi_key\":\"HGRHWEAUHXYNNP-UHFFFAOYSA-N\"},\"molecular_properties\":{\"molecular_weight\":270.376,\"exact_mass\":270.17321332,\"monoisotopic_mass\":270.17321332,\"xlogp3\":3.393,\"clogp\":3.393,\"hbd_count\":1,\"hba_count\":2,\"rotatable_bond_count\":8,\"heavy_atom_count\":20,\"aromatic_ring_count\":2,\"aliphatic_ring_count\":0,\"tpsa\":28.26,\"formal_charge\":0,\"bertz_complexity\":575.688,\"defined_atom_stereocenter_count\":0,\"undefined_atom_stereocenter_count\":0,\"covalently_bonded_unit_count\":1,\"compound_is_canonicalized\":true,\"has_stereochemistry\":false,\"radical_electron_count\":0,\"isotope_atom_count\":0,\"total_charge\":0},\"drug_properties\":{\"lipinski_violations\":0,\"lipinski_mw_violation\":false,\"lipinski_logp_violation\":false,\"lipinski_hbd_violation\":false,\"lipinski_hba_violation\":false,\"veber_violations\":0,\"veber_rotatable_bonds_violation\":false,\"veber_tpsa_violation\":false,\"qed_score\":0.7447,\"drug_like\":true,\"oral_bioavailability_likely\":false},\"physical_properties\":{\"estimated_logS\":-3.9,\"estimated_melting_point\":null,\"estimated_boiling_point\":null},\"additional_descriptors\":{\"kappa1\":14.41,\"kappa2\":6.9632,\"kappa3\":3.3058,\"chi0v\":12.0622,\"chi1v\":6.863,\"chi2v\":4.8124,\"ring_count\":2,\"saturated_ring_count\":0,\"heterocycle_count\":1,\"carbon_count\":17,\"nitrogen_count\":2,\"oxygen_count\":1,\"sulfur_count\":0,\"phosphorus_count\":0,\"halogen_count\":0},\"analysis_metadata\":{\"rdkit_version\":\"2025.03.3\",\"canonical_smiles\":\"C=CCN(CC=C)CCc1c[nH]c2ccc(OC)cc12\",\"analysis_success\":true},\"success\":true}', '2025-07-11 20:27:59', 'rdkit', '2025.03.3');

-- --------------------------------------------------------

--
-- Table structure for table `chemicals`
--

CREATE TABLE `chemicals` (
  `id` int(11) NOT NULL,
  `smiles` varchar(1000) NOT NULL,
  `smiles_canonical` varchar(1000) DEFAULT NULL,
  `iupac_name` text DEFAULT NULL,
  `common_name` varchar(500) DEFAULT NULL,
  `inchi` text DEFAULT NULL,
  `inchi_key` varchar(255) DEFAULT NULL,
  `molecular_formula` varchar(255) DEFAULT NULL,
  `cas_number` varchar(50) DEFAULT NULL,
  `pubchem_cid` bigint(20) DEFAULT NULL,
  `chemspider_id` bigint(20) DEFAULT NULL,
  `created_at` timestamp NOT NULL DEFAULT current_timestamp(),
  `updated_at` timestamp NOT NULL DEFAULT current_timestamp() ON UPDATE current_timestamp()
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

--
-- Dumping data for table `chemicals`
--

INSERT INTO `chemicals` (`id`, `smiles`, `smiles_canonical`, `iupac_name`, `common_name`, `inchi`, `inchi_key`, `molecular_formula`, `cas_number`, `pubchem_cid`, `chemspider_id`, `created_at`, `updated_at`) VALUES
(6, 'CCO', 'CCO', NULL, NULL, 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N', 'C2H6O', NULL, NULL, NULL, '2025-07-11 20:24:12', '2025-07-11 20:24:12'),
(7, 'O', 'O', NULL, NULL, 'InChI=1S/H2O/h1H2', 'XLYOFNOQVPJJNP-UHFFFAOYSA-N', 'H2O', NULL, NULL, NULL, '2025-07-11 20:25:19', '2025-07-11 20:25:19'),
(8, 'C=CCN(CCC1=CNC2=C1C=C(OC)C=C2)CC=C', 'C=CCN(CC=C)CCc1c[nH]c2ccc(OC)cc12', NULL, NULL, 'InChI=1S/C17H22N2O/c1-4-9-19(10-5-2)11-8-14-13-18-17-7-6-15(20-3)12-16(14)17/h4-7,12-13,18H,1-2,8-11H2,3H3', 'HGRHWEAUHXYNNP-UHFFFAOYSA-N', 'C17H22N2O', NULL, NULL, NULL, '2025-07-11 20:27:59', '2025-07-11 20:27:59');

-- --------------------------------------------------------

--
-- Stand-in structure for view `chemical_full_profile`
-- (See below for the actual view)
--
CREATE TABLE `chemical_full_profile` (
`id` int(11)
,`smiles` varchar(1000)
,`smiles_canonical` varchar(1000)
,`iupac_name` text
,`common_name` varchar(500)
,`inchi` text
,`inchi_key` varchar(255)
,`molecular_formula` varchar(255)
,`cas_number` varchar(50)
,`pubchem_cid` bigint(20)
,`chemspider_id` bigint(20)
,`created_at` timestamp
,`updated_at` timestamp
,`molecular_weight` decimal(10,4)
,`exact_mass` decimal(15,8)
,`xlogp3` decimal(8,4)
,`hbd_count` int(11)
,`hba_count` int(11)
,`rotatable_bond_count` int(11)
,`tpsa` decimal(8,3)
,`heavy_atom_count` int(11)
,`formal_charge` int(11)
,`complexity_score` decimal(10,3)
,`lipinski_violations` int(11)
,`sa_score` decimal(4,2)
,`water_solubility` decimal(15,8)
,`melting_point` decimal(8,2)
,`boiling_point` decimal(8,2)
);

-- --------------------------------------------------------

--
-- Table structure for table `chemical_synonyms`
--

CREATE TABLE `chemical_synonyms` (
  `id` int(11) NOT NULL,
  `chemical_id` int(11) NOT NULL,
  `synonym` varchar(500) NOT NULL,
  `synonym_type` varchar(50) DEFAULT NULL,
  `source` varchar(100) DEFAULT NULL,
  `created_at` timestamp NOT NULL DEFAULT current_timestamp()
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `drug_properties`
--

CREATE TABLE `drug_properties` (
  `id` int(11) NOT NULL,
  `chemical_id` int(11) NOT NULL,
  `lipinski_violations` int(11) DEFAULT 0,
  `lipinski_mw_violation` tinyint(1) DEFAULT 0,
  `lipinski_logp_violation` tinyint(1) DEFAULT 0,
  `lipinski_hbd_violation` tinyint(1) DEFAULT 0,
  `lipinski_hba_violation` tinyint(1) DEFAULT 0,
  `veber_violations` int(11) DEFAULT 0,
  `veber_rotatable_bonds_violation` tinyint(1) DEFAULT 0,
  `veber_tpsa_violation` tinyint(1) DEFAULT 0,
  `ghose_violations` int(11) DEFAULT 0,
  `pains_alerts` int(11) DEFAULT 0,
  `pains_description` text DEFAULT NULL,
  `sa_score` decimal(4,2) DEFAULT NULL,
  `qed_score` decimal(4,3) DEFAULT NULL,
  `drug_like` tinyint(1) DEFAULT 0,
  `oral_bioavailability_likely` tinyint(1) DEFAULT 0,
  `bbb_penetration_probability` decimal(4,3) DEFAULT NULL,
  `hia_probability` decimal(4,3) DEFAULT NULL,
  `cyp3a4_substrate_probability` decimal(4,3) DEFAULT NULL,
  `cyp3a4_inhibitor_probability` decimal(4,3) DEFAULT NULL,
  `cyp2d6_substrate_probability` decimal(4,3) DEFAULT NULL,
  `cyp2d6_inhibitor_probability` decimal(4,3) DEFAULT NULL,
  `created_at` timestamp NOT NULL DEFAULT current_timestamp(),
  `updated_at` timestamp NOT NULL DEFAULT current_timestamp() ON UPDATE current_timestamp()
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

--
-- Dumping data for table `drug_properties`
--

INSERT INTO `drug_properties` (`id`, `chemical_id`, `lipinski_violations`, `lipinski_mw_violation`, `lipinski_logp_violation`, `lipinski_hbd_violation`, `lipinski_hba_violation`, `veber_violations`, `veber_rotatable_bonds_violation`, `veber_tpsa_violation`, `ghose_violations`, `pains_alerts`, `pains_description`, `sa_score`, `qed_score`, `drug_like`, `oral_bioavailability_likely`, `bbb_penetration_probability`, `hia_probability`, `cyp3a4_substrate_probability`, `cyp3a4_inhibitor_probability`, `cyp2d6_substrate_probability`, `cyp2d6_inhibitor_probability`, `created_at`, `updated_at`) VALUES
(2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NULL, NULL, 0.407, 1, 1, NULL, NULL, NULL, NULL, NULL, NULL, '2025-07-11 20:24:12', '2025-07-11 20:24:12'),
(3, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NULL, NULL, 0.328, 1, 1, NULL, NULL, NULL, NULL, NULL, NULL, '2025-07-11 20:25:19', '2025-07-11 20:25:19'),
(4, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NULL, NULL, 0.745, 1, 0, NULL, NULL, NULL, NULL, NULL, NULL, '2025-07-11 20:27:59', '2025-07-11 20:27:59');

-- --------------------------------------------------------

--
-- Table structure for table `molecular_properties`
--

CREATE TABLE `molecular_properties` (
  `id` int(11) NOT NULL,
  `chemical_id` int(11) NOT NULL,
  `molecular_weight` decimal(10,4) DEFAULT NULL,
  `exact_mass` decimal(15,8) DEFAULT NULL,
  `monoisotopic_mass` decimal(15,8) DEFAULT NULL,
  `xlogp3` decimal(8,4) DEFAULT NULL,
  `clogp` decimal(8,4) DEFAULT NULL,
  `hbd_count` int(11) DEFAULT NULL,
  `hba_count` int(11) DEFAULT NULL,
  `rotatable_bond_count` int(11) DEFAULT NULL,
  `heavy_atom_count` int(11) DEFAULT NULL,
  `aromatic_ring_count` int(11) DEFAULT NULL,
  `aliphatic_ring_count` int(11) DEFAULT NULL,
  `tpsa` decimal(8,3) DEFAULT NULL,
  `molecular_volume` decimal(10,3) DEFAULT NULL,
  `formal_charge` int(11) DEFAULT NULL,
  `total_charge` int(11) DEFAULT NULL,
  `complexity_score` decimal(10,3) DEFAULT NULL,
  `bertz_complexity` decimal(10,3) DEFAULT NULL,
  `isotope_atom_count` int(11) DEFAULT NULL,
  `radical_electron_count` int(11) DEFAULT NULL,
  `defined_atom_stereocenter_count` int(11) DEFAULT NULL,
  `undefined_atom_stereocenter_count` int(11) DEFAULT NULL,
  `defined_bond_stereocenter_count` int(11) DEFAULT NULL,
  `undefined_bond_stereocenter_count` int(11) DEFAULT NULL,
  `covalently_bonded_unit_count` int(11) DEFAULT NULL,
  `component_count` int(11) DEFAULT NULL,
  `compound_is_canonicalized` tinyint(1) DEFAULT 0,
  `has_stereochemistry` tinyint(1) DEFAULT 0,
  `created_at` timestamp NOT NULL DEFAULT current_timestamp(),
  `updated_at` timestamp NOT NULL DEFAULT current_timestamp() ON UPDATE current_timestamp()
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

--
-- Dumping data for table `molecular_properties`
--

INSERT INTO `molecular_properties` (`id`, `chemical_id`, `molecular_weight`, `exact_mass`, `monoisotopic_mass`, `xlogp3`, `clogp`, `hbd_count`, `hba_count`, `rotatable_bond_count`, `heavy_atom_count`, `aromatic_ring_count`, `aliphatic_ring_count`, `tpsa`, `molecular_volume`, `formal_charge`, `total_charge`, `complexity_score`, `bertz_complexity`, `isotope_atom_count`, `radical_electron_count`, `defined_atom_stereocenter_count`, `undefined_atom_stereocenter_count`, `defined_bond_stereocenter_count`, `undefined_bond_stereocenter_count`, `covalently_bonded_unit_count`, `component_count`, `compound_is_canonicalized`, `has_stereochemistry`, `created_at`, `updated_at`) VALUES
(6, 6, 46.0690, 46.04186481, 46.04186481, -0.0014, -0.0014, 1, 1, 0, 3, 0, 0, 20.230, NULL, 0, 0, NULL, 2.755, 0, 0, 0, 0, NULL, NULL, 1, NULL, 1, 0, '2025-07-11 20:24:12', '2025-07-11 20:24:12'),
(7, 7, 18.0150, 18.01056468, 18.01056468, -0.8247, -0.8247, 0, 0, 0, 1, 0, 0, 31.500, NULL, 0, 0, NULL, 0.000, 0, 0, 0, 0, NULL, NULL, 1, NULL, 1, 0, '2025-07-11 20:25:19', '2025-07-11 20:25:19'),
(8, 8, 270.3760, 270.17321332, 270.17321332, 3.3930, 3.3930, 1, 2, 8, 20, 2, 0, 28.260, NULL, 0, 0, NULL, 575.688, 0, 0, 0, 0, NULL, NULL, 1, NULL, 1, 0, '2025-07-11 20:27:59', '2025-07-11 20:27:59');

-- --------------------------------------------------------

--
-- Table structure for table `physical_properties`
--

CREATE TABLE `physical_properties` (
  `id` int(11) NOT NULL,
  `chemical_id` int(11) NOT NULL,
  `physical_state` varchar(50) DEFAULT NULL,
  `color` varchar(100) DEFAULT NULL,
  `odor` varchar(200) DEFAULT NULL,
  `melting_point` decimal(8,2) DEFAULT NULL,
  `boiling_point` decimal(8,2) DEFAULT NULL,
  `flash_point` decimal(8,2) DEFAULT NULL,
  `autoignition_temperature` decimal(8,2) DEFAULT NULL,
  `density` decimal(8,4) DEFAULT NULL,
  `vapor_pressure` decimal(15,8) DEFAULT NULL,
  `vapor_density` decimal(8,3) DEFAULT NULL,
  `water_solubility` decimal(15,8) DEFAULT NULL,
  `water_solubility_unit` varchar(20) DEFAULT 'mg/L',
  `logS` decimal(6,3) DEFAULT NULL,
  `refractive_index` decimal(8,5) DEFAULT NULL,
  `stability` varchar(500) DEFAULT NULL,
  `shelf_life` varchar(200) DEFAULT NULL,
  `created_at` timestamp NOT NULL DEFAULT current_timestamp(),
  `updated_at` timestamp NOT NULL DEFAULT current_timestamp() ON UPDATE current_timestamp(),
  `estimated_logS` decimal(6,3) DEFAULT NULL,
  `estimated_melting_point` decimal(8,2) DEFAULT NULL,
  `estimated_boiling_point` decimal(8,2) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

--
-- Dumping data for table `physical_properties`
--

INSERT INTO `physical_properties` (`id`, `chemical_id`, `physical_state`, `color`, `odor`, `melting_point`, `boiling_point`, `flash_point`, `autoignition_temperature`, `density`, `vapor_pressure`, `vapor_density`, `water_solubility`, `water_solubility_unit`, `logS`, `refractive_index`, `stability`, `shelf_life`, `created_at`, `updated_at`, `estimated_logS`, `estimated_melting_point`, `estimated_boiling_point`) VALUES
(1, 6, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 'mg/L', NULL, NULL, NULL, NULL, '2025-07-11 20:24:12', '2025-07-11 20:24:12', 0.040, NULL, NULL),
(2, 7, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 'mg/L', NULL, NULL, NULL, NULL, '2025-07-11 20:25:19', '2025-07-11 20:25:19', 0.732, NULL, NULL),
(3, 8, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 'mg/L', NULL, NULL, NULL, NULL, '2025-07-11 20:27:59', '2025-07-11 20:27:59', -3.900, NULL, NULL);

-- --------------------------------------------------------

--
-- Table structure for table `safety_properties`
--

CREATE TABLE `safety_properties` (
  `id` int(11) NOT NULL,
  `chemical_id` int(11) NOT NULL,
  `ghs_pictograms` longtext CHARACTER SET utf8mb4 COLLATE utf8mb4_bin DEFAULT NULL CHECK (json_valid(`ghs_pictograms`)),
  `ghs_signal_word` varchar(50) DEFAULT NULL,
  `ghs_hazard_statements` longtext CHARACTER SET utf8mb4 COLLATE utf8mb4_bin DEFAULT NULL CHECK (json_valid(`ghs_hazard_statements`)),
  `ghs_precautionary_statements` longtext CHARACTER SET utf8mb4 COLLATE utf8mb4_bin DEFAULT NULL CHECK (json_valid(`ghs_precautionary_statements`)),
  `oral_ld50_rat` decimal(12,3) DEFAULT NULL,
  `dermal_ld50_rabbit` decimal(12,3) DEFAULT NULL,
  `inhalation_lc50_rat` decimal(12,3) DEFAULT NULL,
  `ames_test_result` varchar(20) DEFAULT NULL,
  `carcinogenicity_classification` varchar(100) DEFAULT NULL,
  `bioconcentration_factor` decimal(12,3) DEFAULT NULL,
  `biodegradation_probability` decimal(4,3) DEFAULT NULL,
  `incompatible_materials` text DEFAULT NULL,
  `hazardous_decomposition` text DEFAULT NULL,
  `created_at` timestamp NOT NULL DEFAULT current_timestamp(),
  `updated_at` timestamp NOT NULL DEFAULT current_timestamp() ON UPDATE current_timestamp()
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `search_history`
--

CREATE TABLE `search_history` (
  `id` int(11) NOT NULL,
  `session_id` varchar(255) DEFAULT NULL,
  `chemical_id` int(11) DEFAULT NULL,
  `search_query` varchar(1000) DEFAULT NULL,
  `search_type` varchar(50) DEFAULT NULL,
  `ip_address` varchar(45) DEFAULT NULL,
  `user_agent` text DEFAULT NULL,
  `search_date` timestamp NOT NULL DEFAULT current_timestamp()
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

--
-- Dumping data for table `search_history`
--

INSERT INTO `search_history` (`id`, `session_id`, `chemical_id`, `search_query`, `search_type`, `ip_address`, `user_agent`, `search_date`) VALUES
(1, 'chm02q2sgusrms869nfkoigu03', 6, 'CCO', 'smiles', '127.0.0.1', 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:140.0) Gecko/20100101 Firefox/140.0', '2025-07-11 20:24:12'),
(2, 'n3puqcqepdb8a0nu4pomqchhr1', 7, 'O', 'smiles', '::1', 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Code-Insiders/1.103.0-insider Chrome/134.0.6998.205 Electron/35.6.0 Safari/537.36', '2025-07-11 20:25:20'),
(3, 'chm02q2sgusrms869nfkoigu03', 8, 'C=CCN(CCC1=CNC2=C1C=C(OC)C=C2)CC=C', 'smiles', '127.0.0.1', 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:140.0) Gecko/20100101 Firefox/140.0', '2025-07-11 20:27:59');

-- --------------------------------------------------------

--
-- Structure for view `chemical_full_profile`
--
DROP TABLE IF EXISTS `chemical_full_profile`;

CREATE ALGORITHM=UNDEFINED DEFINER=`root`@`localhost` SQL SECURITY DEFINER VIEW `chemical_full_profile`  AS SELECT `c`.`id` AS `id`, `c`.`smiles` AS `smiles`, `c`.`smiles_canonical` AS `smiles_canonical`, `c`.`iupac_name` AS `iupac_name`, `c`.`common_name` AS `common_name`, `c`.`inchi` AS `inchi`, `c`.`inchi_key` AS `inchi_key`, `c`.`molecular_formula` AS `molecular_formula`, `c`.`cas_number` AS `cas_number`, `c`.`pubchem_cid` AS `pubchem_cid`, `c`.`chemspider_id` AS `chemspider_id`, `c`.`created_at` AS `created_at`, `c`.`updated_at` AS `updated_at`, `mp`.`molecular_weight` AS `molecular_weight`, `mp`.`exact_mass` AS `exact_mass`, `mp`.`xlogp3` AS `xlogp3`, `mp`.`hbd_count` AS `hbd_count`, `mp`.`hba_count` AS `hba_count`, `mp`.`rotatable_bond_count` AS `rotatable_bond_count`, `mp`.`tpsa` AS `tpsa`, `mp`.`heavy_atom_count` AS `heavy_atom_count`, `mp`.`formal_charge` AS `formal_charge`, `mp`.`complexity_score` AS `complexity_score`, `dp`.`lipinski_violations` AS `lipinski_violations`, `dp`.`sa_score` AS `sa_score`, `pp`.`water_solubility` AS `water_solubility`, `pp`.`melting_point` AS `melting_point`, `pp`.`boiling_point` AS `boiling_point` FROM (((`chemicals` `c` left join `molecular_properties` `mp` on(`c`.`id` = `mp`.`chemical_id`)) left join `drug_properties` `dp` on(`c`.`id` = `dp`.`chemical_id`)) left join `physical_properties` `pp` on(`c`.`id` = `pp`.`chemical_id`)) ;

--
-- Indexes for dumped tables
--

--
-- Indexes for table `analysis_cache`
--
ALTER TABLE `analysis_cache`
  ADD PRIMARY KEY (`id`),
  ADD KEY `chemical_id` (`chemical_id`),
  ADD KEY `idx_analysis_type` (`analysis_type`),
  ADD KEY `idx_analysis_date` (`analysis_date`);

--
-- Indexes for table `chemicals`
--
ALTER TABLE `chemicals`
  ADD PRIMARY KEY (`id`),
  ADD UNIQUE KEY `idx_smiles_canonical` (`smiles_canonical`) USING HASH,
  ADD KEY `idx_inchi_key` (`inchi_key`),
  ADD KEY `idx_cas_number` (`cas_number`),
  ADD KEY `idx_pubchem_cid` (`pubchem_cid`),
  ADD KEY `idx_molecular_formula` (`molecular_formula`);

--
-- Indexes for table `chemical_synonyms`
--
ALTER TABLE `chemical_synonyms`
  ADD PRIMARY KEY (`id`),
  ADD KEY `chemical_id` (`chemical_id`),
  ADD KEY `idx_synonym` (`synonym`),
  ADD KEY `idx_synonym_type` (`synonym_type`);

--
-- Indexes for table `drug_properties`
--
ALTER TABLE `drug_properties`
  ADD PRIMARY KEY (`id`),
  ADD KEY `chemical_id` (`chemical_id`);

--
-- Indexes for table `molecular_properties`
--
ALTER TABLE `molecular_properties`
  ADD PRIMARY KEY (`id`),
  ADD KEY `chemical_id` (`chemical_id`),
  ADD KEY `idx_molecular_weight` (`molecular_weight`),
  ADD KEY `idx_xlogp3` (`xlogp3`),
  ADD KEY `idx_tpsa` (`tpsa`);

--
-- Indexes for table `physical_properties`
--
ALTER TABLE `physical_properties`
  ADD PRIMARY KEY (`id`),
  ADD KEY `chemical_id` (`chemical_id`);

--
-- Indexes for table `safety_properties`
--
ALTER TABLE `safety_properties`
  ADD PRIMARY KEY (`id`),
  ADD KEY `chemical_id` (`chemical_id`);

--
-- Indexes for table `search_history`
--
ALTER TABLE `search_history`
  ADD PRIMARY KEY (`id`),
  ADD KEY `chemical_id` (`chemical_id`),
  ADD KEY `idx_session_id` (`session_id`),
  ADD KEY `idx_search_date` (`search_date`),
  ADD KEY `idx_search_type` (`search_type`);

--
-- AUTO_INCREMENT for dumped tables
--

--
-- AUTO_INCREMENT for table `analysis_cache`
--
ALTER TABLE `analysis_cache`
  MODIFY `id` int(11) NOT NULL AUTO_INCREMENT, AUTO_INCREMENT=4;

--
-- AUTO_INCREMENT for table `chemicals`
--
ALTER TABLE `chemicals`
  MODIFY `id` int(11) NOT NULL AUTO_INCREMENT, AUTO_INCREMENT=9;

--
-- AUTO_INCREMENT for table `chemical_synonyms`
--
ALTER TABLE `chemical_synonyms`
  MODIFY `id` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `drug_properties`
--
ALTER TABLE `drug_properties`
  MODIFY `id` int(11) NOT NULL AUTO_INCREMENT, AUTO_INCREMENT=5;

--
-- AUTO_INCREMENT for table `molecular_properties`
--
ALTER TABLE `molecular_properties`
  MODIFY `id` int(11) NOT NULL AUTO_INCREMENT, AUTO_INCREMENT=9;

--
-- AUTO_INCREMENT for table `physical_properties`
--
ALTER TABLE `physical_properties`
  MODIFY `id` int(11) NOT NULL AUTO_INCREMENT, AUTO_INCREMENT=4;

--
-- AUTO_INCREMENT for table `safety_properties`
--
ALTER TABLE `safety_properties`
  MODIFY `id` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT for table `search_history`
--
ALTER TABLE `search_history`
  MODIFY `id` int(11) NOT NULL AUTO_INCREMENT, AUTO_INCREMENT=4;

--
-- Constraints for dumped tables
--

--
-- Constraints for table `analysis_cache`
--
ALTER TABLE `analysis_cache`
  ADD CONSTRAINT `analysis_cache_ibfk_1` FOREIGN KEY (`chemical_id`) REFERENCES `chemicals` (`id`) ON DELETE CASCADE;

--
-- Constraints for table `chemical_synonyms`
--
ALTER TABLE `chemical_synonyms`
  ADD CONSTRAINT `chemical_synonyms_ibfk_1` FOREIGN KEY (`chemical_id`) REFERENCES `chemicals` (`id`) ON DELETE CASCADE;

--
-- Constraints for table `drug_properties`
--
ALTER TABLE `drug_properties`
  ADD CONSTRAINT `drug_properties_ibfk_1` FOREIGN KEY (`chemical_id`) REFERENCES `chemicals` (`id`) ON DELETE CASCADE;

--
-- Constraints for table `molecular_properties`
--
ALTER TABLE `molecular_properties`
  ADD CONSTRAINT `molecular_properties_ibfk_1` FOREIGN KEY (`chemical_id`) REFERENCES `chemicals` (`id`) ON DELETE CASCADE;

--
-- Constraints for table `physical_properties`
--
ALTER TABLE `physical_properties`
  ADD CONSTRAINT `physical_properties_ibfk_1` FOREIGN KEY (`chemical_id`) REFERENCES `chemicals` (`id`) ON DELETE CASCADE;

--
-- Constraints for table `safety_properties`
--
ALTER TABLE `safety_properties`
  ADD CONSTRAINT `safety_properties_ibfk_1` FOREIGN KEY (`chemical_id`) REFERENCES `chemicals` (`id`) ON DELETE CASCADE;

--
-- Constraints for table `search_history`
--
ALTER TABLE `search_history`
  ADD CONSTRAINT `search_history_ibfk_1` FOREIGN KEY (`chemical_id`) REFERENCES `chemicals` (`id`) ON DELETE SET NULL;
COMMIT;

/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
