/**
 * Copyright (c) 2017-present, Facebook, Inc.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

const React = require('react');

const CompLibrary = require('../../core/CompLibrary.js');

const Container = CompLibrary.Container;
const GridBlock = CompLibrary.GridBlock;

function Help(props) {
  const { config: siteConfig, language = '' } = props;
  const { baseUrl, docsUrl } = siteConfig;
  const docsPart = `${docsUrl ? `${docsUrl}/` : ''}`;
  const langPart = `${language ? `${language}/` : ''}`;
  const docUrl = doc => `${baseUrl}${docsPart}${langPart}${doc}`;

  const supportLinks = [
    {
      content: `Learn more using the [documentation on this site.](${docUrl(
        'intro_public.html',
      )})`,
      title: 'Browse Docs',
    },
    {
      content: 'Feel free to ask questions or raise issues on the [GitHub Repository](https://www.github.com/julianstanley/ProteinFeatures).',
      title: 'Interact via GitHub',
    },
    {
      content: "For inquiries, feel free to contact [Julian](mailto:julianstanleya@gmail.com), who created this documentation, or [Dr. Roger Chang](roger_chang@hms.harvard.edu), who spearheaded the project.",
      title: 'Contact us Directly',
    },
  ];

  return (
    <div className="docMainWrapper wrapper">
      <Container className="mainContainer documentContainer postContainer">
        <div className="post">
          <header className="postHeader">
            <h1>Need help?</h1>
          </header>
          <p>This project is maintained by the Silver Lab at Harvard Medical School.</p>
          <p>
          </p>
          <GridBlock contents={supportLinks} layout="threeColumn" />
        </div>
      </Container>
    </div>
  );
}

module.exports = Help;
