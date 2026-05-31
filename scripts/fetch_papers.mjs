#!/usr/bin/env node
// Fetch reference papers via headed Playwright browser (institutional/TIB access required).
// Usage: node scripts/fetch_papers.mjs
//
// 1. Ensure institutional access is active (for example, TIB VPN).
// 2. Script opens headed Chromium with persistent profile.
// 3. Solve any Cloudflare challenges in the browser window when prompted.
// 4. Script fetches PDFs via page.request.get() after challenge is cleared.

import { chromium } from 'playwright';
import { writeFileSync, existsSync, mkdirSync } from 'fs';
import { resolve } from 'path';

const PAPERS_DIR = resolve(import.meta.dirname, '..', 'docs', 'papers');
const PROFILE_DIR = resolve(PAPERS_DIR, '.browser-profile');

// Papers grouped by publisher domain for efficient Cloudflare handling.
// Direct PDF URLs where known; DOI fallback otherwise.
const PAPERS = [
  // OUP (Oxford University Press)
  {
    label: 'DecentKing2008',
    desc: 'Decent & King 2008 — main paper, slender cone recoil, IMAJAM 73(1), 37-68',
    pdfUrl: null,
    doiFallback: '10.1093/imamat/hxm043',
    file: 'DecentKing2008_IMAJAM_73_1_37-68_hxm043.pdf',
    domain: 'oup',
  },
  // SIAM
  {
    label: 'KellerMiksis1983',
    desc: 'Keller & Miksis 1983 — t^{2/3} similarity scaling',
    abstractUrl: 'https://epubs.siam.org/doi/10.1137/0143018',
    pdfUrl: 'https://epubs.siam.org/doi/pdf/10.1137/0143018',
    file: 'KellerMiksis1983_SIAMJAM_43_268.pdf',
    domain: 'siam',
  },
  // Cambridge University Press (JFM)
  {
    label: 'Billingham1999',
    desc: 'Billingham 1999 — fat cone with Hankel-Laplace',
    abstractUrl: 'https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/surfacetensiondriven-flow-in-fat-fluid-wedges-and-cones/47EE570853A6F2CCA2848E3B17C8B5A8',
    pdfUrl: null, // will resolve via DOI
    doiFallback: '10.1017/S0022112099006047',
    file: 'Billingham1999_JFM_397_45.pdf',
    domain: 'cambridge',
  },
  // APS (American Physical Society) — already downloaded
  {
    label: 'Eggers1997',
    desc: 'Eggers 1997 — review with context',
    pdfUrl: 'https://journals.aps.org/rmp/pdf/10.1103/RevModPhys.69.865',
    file: 'Eggers1997_RMP_69_865.pdf',
    domain: 'aps',
  },
  // AIP (Physics of Fluids)
  {
    label: 'KellerKingTing1995',
    desc: 'Keller, King & Ting 1995 — blob formation',
    abstractUrl: 'https://pubs.aip.org/aip/pof/article-abstract/7/1/226/262287',
    pdfUrl: null, // will resolve via DOI
    doiFallback: '10.1063/1.868723',
    file: 'KellerKingTing1995_PoF_7_226.pdf',
    domain: 'aip',
  },
  // Springer (IUTAM proceedings)
  {
    label: 'DecentKing2001',
    desc: 'Decent & King 2001 — conference proceedings, recoil of broken liquid bridge',
    pdfUrl: null,
    doiFallback: '10.1007/978-94-010-0796-2_10',
    file: 'DecentKing2001_IUTAM.pdf',
    domain: 'springer',
  },
];

async function waitForChallenge(page, label, timeoutMs = 180000) {
  // Check if we're on a Cloudflare challenge page
  const content = await page.content();
  const isChallenge = content.includes('challenge-platform')
    || content.includes('cf-challenge')
    || content.includes('Just a moment')
    || content.includes('Checking your browser');

  if (!isChallenge) return true;

  console.log(`\n>>> Cloudflare challenge detected for ${label} <<<`);
  console.log(`>>> Click the challenge in the browser window <<<`);
  console.log(`>>> You have ${timeoutMs / 1000}s <<<\n`);

  try {
    await page.waitForFunction(() => {
      return !document.title.includes('Just a moment')
        && !document.querySelector('#challenge-running')
        && !document.querySelector('.cf-challenge');
    }, { timeout: timeoutMs });
    console.log('Challenge passed!');
    await new Promise(r => setTimeout(r, 3000));
    return true;
  } catch {
    console.log('Challenge timed out.');
    return false;
  }
}

async function fetchPDF(page, paper) {
  const outPath = resolve(PAPERS_DIR, paper.file);
  if (existsSync(outPath)) {
    console.log(`SKIP ${paper.label}: ${paper.file} (already exists)`);
    return 'skip';
  }

  console.log(`\n=== ${paper.label}: ${paper.desc} ===`);

  // Strategy 1: Direct PDF URL if known
  if (paper.pdfUrl) {
    const result = await tryDirectDownload(page, paper.pdfUrl, outPath, paper.label);
    if (result) return 'ok';
  }

  // Strategy 2: Navigate to abstract page, find PDF link
  if (paper.abstractUrl) {
    const result = await tryViaAbstract(page, paper.abstractUrl, outPath, paper.label);
    if (result) return 'ok';
  }

  // Strategy 3: Resolve DOI, find PDF on publisher page
  if (paper.doiFallback) {
    const result = await tryViaDOI(page, paper.doiFallback, outPath, paper.label);
    if (result) return 'ok';
  }

  console.log(`FAIL ${paper.label}: all strategies exhausted`);
  return 'fail';
}

async function tryDirectDownload(page, url, outPath, label) {
  try {
    process.stdout.write(`  Download: ${url} ... `);

    // Navigate the browser itself to the PDF URL so institutional cookies are sent.
    // Intercept the response to capture the PDF body.
    const responsePromise = page.waitForResponse(
      resp => resp.url().includes(new URL(url).pathname) || resp.url() === url,
      { timeout: 30000 }
    ).catch(() => null);

    const response = await page.goto(url, { waitUntil: 'commit', timeout: 30000 }).catch(() => null);
    const resp = response || await responsePromise;

    if (!resp) {
      console.log('No response');
      return false;
    }

    if (resp.status() !== 200) {
      console.log(`HTTP ${resp.status()}`);
      return false;
    }

    const body = await resp.body();
    const header = body.slice(0, 5).toString();
    if (header !== '%PDF-') {
      console.log(`Not PDF (got: "${header.replace(/\n/g, '\\n')}")`);
      // Maybe we landed on an HTML page with a PDF link
      return false;
    }

    writeFileSync(outPath, body);
    console.log(`OK (${(body.length / 1024).toFixed(0)} KB)`);
    return true;
  } catch (e) {
    console.log(`Error: ${e.message}`);
    return false;
  }
}

async function tryViaAbstract(page, url, outPath, label) {
  try {
    console.log(`  Abstract: ${url}`);
    await page.goto(url, { waitUntil: 'domcontentloaded', timeout: 30000 });
    await waitForChallenge(page, label);

    // Look for PDF links on the abstract page
    const pdfHref = await findPDFHref(page);
    if (!pdfHref) {
      console.log(`  No PDF link found on abstract page.`);
      return false;
    }

    const fullUrl = pdfHref.startsWith('http') ? pdfHref : new URL(pdfHref, page.url()).href;
    return await tryDirectDownload(page, fullUrl, outPath, label);
  } catch (e) {
    console.log(`  Abstract error: ${e.message}`);
    return false;
  }
}

async function tryViaDOI(page, doi, outPath, label) {
  try {
    const doiUrl = `https://doi.org/${doi}`;
    console.log(`  DOI: ${doiUrl}`);
    await page.goto(doiUrl, { waitUntil: 'domcontentloaded', timeout: 30000 });
    await waitForChallenge(page, label);

    const landed = page.url();
    console.log(`  Landed: ${landed}`);

    // Try PDF link on the resolved page
    const pdfHref = await findPDFHref(page);
    if (pdfHref) {
      const fullUrl = pdfHref.startsWith('http') ? pdfHref : new URL(pdfHref, page.url()).href;
      return await tryDirectDownload(page, fullUrl, outPath, label);
    }

    // Publisher-specific PDF URL patterns
    if (landed.includes('cambridge.org')) {
      // Cambridge: /core/.../<id> → /core/services/aop-cambridge-core/content/view/<id>
      const match = landed.match(/\/article\/([^/]+)\/([a-f0-9]+)/i);
      if (match) {
        const pdfTry = landed.replace('/article/', '/article/pdf/');
        return await tryDirectDownload(page, pdfTry, outPath, label);
      }
    }

    if (landed.includes('link.springer.com')) {
      // Springer: chapter page → /content/pdf/DOI.pdf
      const pdfTry = `https://link.springer.com/content/pdf/${doi}.pdf`;
      return await tryDirectDownload(page, pdfTry, outPath, label);
    }

    if (landed.includes('pubs.aip.org')) {
      // AIP: try appending /pdf or using the PDF endpoint
      const pdfTry = landed.replace('/article/', '/article-pdf/');
      return await tryDirectDownload(page, pdfTry, outPath, label);
    }

    console.log(`  No PDF strategy for: ${landed}`);
    return false;
  } catch (e) {
    console.log(`  DOI error: ${e.message}`);
    return false;
  }
}

async function findPDFHref(page) {
  const selectors = [
    'a[href*="/pdf/"]',
    'a[href$=".pdf"]',
    'a:has-text("Download PDF")',
    'a:has-text("Full Text PDF")',
    'a:has-text("Get PDF")',
    'a:has-text("PDF")',
    'a[title*="PDF"]',
    'a[aria-label*="PDF"]',
    'a.al-link.pdf',
    'a[class*="pdf"]',
    'a[data-article-pdf]',
  ];

  for (const sel of selectors) {
    try {
      const el = await page.$(sel);
      if (el) {
        const href = await el.getAttribute('href');
        if (href && !href.includes('login') && !href.includes('signup')) {
          return href;
        }
      }
    } catch {}
  }
  return null;
}

async function main() {
  mkdirSync(PAPERS_DIR, { recursive: true });
  mkdirSync(PROFILE_DIR, { recursive: true });

  console.log('Launching HEADED Chromium (persistent profile)...');
  console.log('Make sure LUH VPN is active!\n');

  const context = await chromium.launchPersistentContext(PROFILE_DIR, {
    headless: false,
    args: ['--disable-blink-features=AutomationControlled'],
    viewport: { width: 900, height: 600 },
  });
  const page = context.pages()[0] || await context.newPage();

  const results = {};

  for (const paper of PAPERS) {
    results[paper.label] = await fetchPDF(page, paper);
    await new Promise(r => setTimeout(r, 2000));
  }

  console.log('\n=== RESULTS ===');
  for (const [label, status] of Object.entries(results)) {
    const icon = status === 'ok' ? '+' : status === 'skip' ? '~' : 'X';
    console.log(`[${icon}] ${label}: ${status}`);
  }

  const failed = Object.entries(results).filter(([, s]) => s === 'fail');
  if (failed.length > 0) {
    console.log(`\n${failed.length} papers failed. Browser stays open 120s for manual download...`);
    console.log('Save PDFs to: ' + PAPERS_DIR);
    await new Promise(r => setTimeout(r, 120000));
  }

  await context.close();
}

main().catch(e => { console.error(e); process.exit(1); });
