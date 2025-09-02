// Pricing Page with Stripe Integration
import React, { useState } from 'react';
import {
  Box,
  Container,
  Grid,
  Card,
  CardContent,
  CardActions,
  Typography,
  Button,
  List,
  ListItem,
  ListItemIcon,
  ListItemText,
  Chip,
  Switch,
  FormControlLabel,
  Paper,
  Divider,
  useTheme,
  ThemeProvider,
  createTheme
} from '@mui/material';
import {
  Check as CheckIcon,
  Close as CloseIcon,
  Star as StarIcon,
  School as SchoolIcon,
  Business as BusinessIcon,
  Science as ScienceIcon,
  Rocket as RocketIcon,
  Support as SupportIcon,
  Storage as StorageIcon,
  Speed as SpeedIcon,
  Group as GroupIcon
} from '@mui/icons-material';
import { motion } from 'framer-motion';
import Head from 'next/head';
import { loadStripe } from '@stripe/stripe-js';

// Initialize Stripe
const stripePromise = loadStripe(process.env.NEXT_PUBLIC_STRIPE_PUBLIC_KEY || '');

interface PricingTier {
  name: string;
  price: number;
  annual_price: number;
  description: string;
  icon: React.ReactNode;
  color: string;
  popular?: boolean;
  features: string[];
  limits: {
    api_calls: string;
    agents: number;
    pipelines: string;
    compounds: string;
    storage: string;
    support: string;
  };
  cta: string;
  stripe_price_id?: string;
}

const pricingTiers: PricingTier[] = [
  {
    name: 'Free',
    price: 0,
    annual_price: 0,
    description: 'Get started with basic drug discovery tools',
    icon: <ScienceIcon />,
    color: '#9e9e9e',
    features: [
      'Access to 10 PhD agents',
      'Basic molecule viewer',
      'Simple similarity search',
      '100 API calls per day',
      'Community support',
      'Public compound libraries'
    ],
    limits: {
      api_calls: '100/day',
      agents: 10,
      pipelines: '1/month',
      compounds: '100/run',
      storage: '1 GB',
      support: 'Community'
    },
    cta: 'Start Free'
  },
  {
    name: 'Researcher',
    price: 99,
    annual_price: 990,
    description: 'Professional tools for individual researchers',
    icon: <StarIcon />,
    color: '#00e5ff',
    popular: true,
    features: [
      'Access to 50 PhD agents',
      'Advanced similarity search',
      'Clustering algorithms',
      'Export to SDF/CSV',
      'Email support',
      'Private compound libraries',
      'Custom descriptors',
      'Batch processing'
    ],
    limits: {
      api_calls: '5,000/day',
      agents: 50,
      pipelines: '10/month',
      compounds: '10,000/run',
      storage: '100 GB',
      support: 'Email'
    },
    cta: 'Start 14-Day Trial',
    stripe_price_id: 'price_researcher_monthly'
  },
  {
    name: 'Professional',
    price: 499,
    annual_price: 4990,
    description: 'Complete platform for drug discovery teams',
    icon: <RocketIcon />,
    color: '#76ff03',
    features: [
      'All 194 PhD agents',
      'Custom pipelines',
      'API access',
      'Team collaboration',
      'Priority support',
      'Advanced analytics',
      'Webhook integrations',
      'Custom models',
      'Parallel processing',
      'Compliance reports'
    ],
    limits: {
      api_calls: '50,000/day',
      agents: 194,
      pipelines: '100/month',
      compounds: '100,000/run',
      storage: '1 TB',
      support: 'Priority'
    },
    cta: 'Start 14-Day Trial',
    stripe_price_id: 'price_professional_monthly'
  },
  {
    name: 'Enterprise',
    price: 2499,
    annual_price: 24990,
    description: 'Unlimited access with dedicated support',
    icon: <BusinessIcon />,
    color: '#ffd740',
    features: [
      'Everything in Professional',
      'Unlimited API calls',
      'Unlimited pipelines',
      'White-label options',
      'On-premise deployment',
      'Custom agent training',
      'Dedicated account manager',
      'SLA guarantee',
      'Custom integrations',
      'Audit logs',
      'SSO/SAML',
      'Advanced security'
    ],
    limits: {
      api_calls: 'Unlimited',
      agents: 194,
      pipelines: 'Unlimited',
      compounds: 'Unlimited',
      storage: 'Unlimited',
      support: 'Dedicated'
    },
    cta: 'Contact Sales',
    stripe_price_id: 'price_enterprise_monthly'
  },
  {
    name: 'Academic',
    price: 49,
    annual_price: 490,
    description: 'Special pricing for academic institutions',
    icon: <SchoolIcon />,
    color: '#b388ff',
    features: [
      'All 194 PhD agents',
      'Academic collaboration tools',
      'Citation exports',
      'Student accounts',
      'Educational resources',
      'Research data sharing',
      'Grant proposal support'
    ],
    limits: {
      api_calls: '10,000/day',
      agents: 194,
      pipelines: '50/month',
      compounds: '50,000/run',
      storage: '500 GB',
      support: 'Email'
    },
    cta: 'Verify Academic Status',
    stripe_price_id: 'price_academic_monthly'
  }
];

const theme = createTheme({
  palette: {
    mode: 'dark',
    primary: { main: '#00e5ff' },
    secondary: { main: '#76ff03' },
    background: {
      default: '#0a0e27',
      paper: '#151934'
    }
  }
});

export default function Pricing() {
  const [billingPeriod, setBillingPeriod] = useState<'monthly' | 'annual'>('monthly');
  
  const handleSubscribe = async (tier: PricingTier) => {
    if (tier.price === 0) {
      // Redirect to signup for free tier
      window.location.href = '/signup';
      return;
    }
    
    if (tier.name === 'Enterprise') {
      // Redirect to contact form for enterprise
      window.location.href = '/contact-sales';
      return;
    }
    
    if (tier.name === 'Academic') {
      // Redirect to academic verification
      window.location.href = '/academic-verification';
      return;
    }
    
    // Create Stripe checkout session
    try {
      const response = await fetch('/api/create-checkout-session', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          price_id: tier.stripe_price_id,
          tier: tier.name.toLowerCase()
        })
      });
      
      const { sessionUrl } = await response.json();
      window.location.href = sessionUrl;
    } catch (error) {
      console.error('Error creating checkout session:', error);
    }
  };

  return (
    <ThemeProvider theme={theme}>
      <Head>
        <title>Pricing - CriOS Dr. Crowe Coder</title>
        <meta name="description" content="Choose the perfect plan for your drug discovery needs" />
      </Head>

      <Container maxWidth="xl" sx={{ py: 8 }}>
        {/* Header */}
        <motion.div
          initial={{ opacity: 0, y: -20 }}
          animate={{ opacity: 1, y: 0 }}
        >
          <Box sx={{ textAlign: 'center', mb: 6 }}>
            <Typography variant="h2" gutterBottom sx={{
              background: 'linear-gradient(45deg, #00e5ff 30%, #76ff03 90%)',
              WebkitBackgroundClip: 'text',
              WebkitTextFillColor: 'transparent',
              fontWeight: 800
            }}>
              Simple, Transparent Pricing
            </Typography>
            <Typography variant="h5" color="textSecondary" sx={{ mb: 4 }}>
              Accelerate drug discovery with 194 PhD AI agents
            </Typography>
            
            {/* Billing Toggle */}
            <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', gap: 2 }}>
              <Typography>Monthly</Typography>
              <Switch
                checked={billingPeriod === 'annual'}
                onChange={(e) => setBillingPeriod(e.target.checked ? 'annual' : 'monthly')}
                color="primary"
              />
              <Typography>
                Annual
                <Chip label="Save 20%" size="small" color="success" sx={{ ml: 1 }} />
              </Typography>
            </Box>
          </Box>
        </motion.div>

        {/* Pricing Cards */}
        <Grid container spacing={3} sx={{ mb: 8 }}>
          {pricingTiers.map((tier, index) => (
            <Grid item xs={12} md={tier.name === 'Academic' ? 12 : 4} lg={tier.name === 'Academic' ? 12 : 2.4} key={tier.name}>
              <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: index * 0.1 }}
                style={{ height: '100%' }}
              >
                <Card sx={{
                  height: '100%',
                  position: 'relative',
                  border: tier.popular ? `2px solid ${tier.color}` : '1px solid rgba(255,255,255,0.1)',
                  background: tier.popular 
                    ? `linear-gradient(135deg, ${tier.color}10 0%, transparent 100%)`
                    : undefined
                }}>
                  {tier.popular && (
                    <Chip
                      label="MOST POPULAR"
                      size="small"
                      sx={{
                        position: 'absolute',
                        top: -12,
                        left: '50%',
                        transform: 'translateX(-50%)',
                        bgcolor: tier.color,
                        color: 'black',
                        fontWeight: 'bold'
                      }}
                    />
                  )}
                  
                  <CardContent sx={{ pb: 0 }}>
                    <Box sx={{ textAlign: 'center', mb: 3 }}>
                      <Box sx={{ color: tier.color, mb: 2 }}>{tier.icon}</Box>
                      <Typography variant="h5" gutterBottom>{tier.name}</Typography>
                      <Typography variant="body2" color="textSecondary" sx={{ mb: 2 }}>
                        {tier.description}
                      </Typography>
                      
                      <Typography variant="h3" sx={{ fontWeight: 'bold' }}>
                        ${billingPeriod === 'annual' 
                          ? Math.floor(tier.annual_price / 12)
                          : tier.price}
                        <Typography variant="body2" component="span" color="textSecondary">
                          /month
                        </Typography>
                      </Typography>
                      
                      {billingPeriod === 'annual' && tier.price > 0 && (
                        <Typography variant="body2" color="success.main">
                          ${tier.annual_price} billed annually
                        </Typography>
                      )}
                    </Box>
                    
                    <Divider sx={{ my: 2 }} />
                    
                    {/* Features */}
                    <List dense>
                      {tier.features.map((feature) => (
                        <ListItem key={feature} sx={{ px: 0 }}>
                          <ListItemIcon sx={{ minWidth: 32 }}>
                            <CheckIcon sx={{ color: tier.color, fontSize: 18 }} />
                          </ListItemIcon>
                          <ListItemText 
                            primary={feature}
                            primaryTypographyProps={{ variant: 'body2' }}
                          />
                        </ListItem>
                      ))}
                    </List>
                    
                    <Divider sx={{ my: 2 }} />
                    
                    {/* Limits */}
                    <Box sx={{ mt: 2 }}>
                      <Typography variant="body2" color="textSecondary" gutterBottom>
                        Limits:
                      </Typography>
                      <Grid container spacing={1}>
                        <Grid item xs={6}>
                          <Typography variant="caption" color="textSecondary">
                            API Calls: {tier.limits.api_calls}
                          </Typography>
                        </Grid>
                        <Grid item xs={6}>
                          <Typography variant="caption" color="textSecondary">
                            Agents: {tier.limits.agents}
                          </Typography>
                        </Grid>
                        <Grid item xs={6}>
                          <Typography variant="caption" color="textSecondary">
                            Pipelines: {tier.limits.pipelines}
                          </Typography>
                        </Grid>
                        <Grid item xs={6}>
                          <Typography variant="caption" color="textSecondary">
                            Storage: {tier.limits.storage}
                          </Typography>
                        </Grid>
                      </Grid>
                    </Box>
                  </CardContent>
                  
                  <CardActions sx={{ p: 2, pt: 0 }}>
                    <Button
                      fullWidth
                      variant={tier.popular ? "contained" : "outlined"}
                      onClick={() => handleSubscribe(tier)}
                      sx={{
                        bgcolor: tier.popular ? tier.color : undefined,
                        color: tier.popular ? 'black' : tier.color,
                        borderColor: tier.color,
                        '&:hover': {
                          bgcolor: tier.popular ? tier.color : `${tier.color}20`,
                          borderColor: tier.color
                        }
                      }}
                    >
                      {tier.cta}
                    </Button>
                  </CardActions>
                </Card>
              </motion.div>
            </Grid>
          ))}
        </Grid>

        {/* Comparison Table */}
        <Paper sx={{ p: 4, mb: 8 }}>
          <Typography variant="h4" gutterBottom sx={{ textAlign: 'center', mb: 4 }}>
            Detailed Feature Comparison
          </Typography>
          
          <Box sx={{ overflowX: 'auto' }}>
            <table style={{ width: '100%', borderCollapse: 'collapse' }}>
              <thead>
                <tr>
                  <th style={{ padding: 16, textAlign: 'left', borderBottom: '2px solid #333' }}>
                    Feature
                  </th>
                  {pricingTiers.filter(t => t.name !== 'Academic').map(tier => (
                    <th key={tier.name} style={{ 
                      padding: 16, 
                      textAlign: 'center', 
                      borderBottom: '2px solid #333',
                      color: tier.color 
                    }}>
                      {tier.name}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {[
                  { feature: 'PhD AI Agents', values: ['10', '50', '194', '194'] },
                  { feature: 'API Calls/Day', values: ['100', '5,000', '50,000', 'Unlimited'] },
                  { feature: 'Discovery Pipelines', values: ['Basic', 'Advanced', 'All 8', 'Custom'] },
                  { feature: 'Compound Library Size', values: ['1K', '100K', '10M', 'Unlimited'] },
                  { feature: 'Real-time Monitoring', values: [false, true, true, true] },
                  { feature: '3D Molecule Viewer', values: [true, true, true, true] },
                  { feature: 'Similarity Search', values: ['Basic', 'Advanced', 'Advanced', 'Custom'] },
                  { feature: 'Clustering Algorithms', values: [false, true, true, true] },
                  { feature: 'Custom Pipelines', values: [false, false, true, true] },
                  { feature: 'API Access', values: [false, false, true, true] },
                  { feature: 'Team Collaboration', values: [false, false, true, true] },
                  { feature: 'White Label', values: [false, false, false, true] },
                  { feature: 'On-Premise Deploy', values: [false, false, false, true] },
                  { feature: 'Support', values: ['Community', 'Email', 'Priority', 'Dedicated'] }
                ].map((row) => (
                  <tr key={row.feature}>
                    <td style={{ padding: 12, borderBottom: '1px solid #222' }}>
                      {row.feature}
                    </td>
                    {row.values.map((value, i) => (
                      <td key={i} style={{ 
                        padding: 12, 
                        textAlign: 'center', 
                        borderBottom: '1px solid #222' 
                      }}>
                        {typeof value === 'boolean' ? (
                          value ? <CheckIcon color="success" /> : <CloseIcon color="error" />
                        ) : (
                          value
                        )}
                      </td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
          </Box>
        </Paper>

        {/* FAQ */}
        <Box sx={{ mb: 8 }}>
          <Typography variant="h4" gutterBottom sx={{ textAlign: 'center', mb: 4 }}>
            Frequently Asked Questions
          </Typography>
          
          <Grid container spacing={3}>
            {[
              {
                q: 'Can I change plans anytime?',
                a: 'Yes! You can upgrade or downgrade your plan at any time. Changes take effect immediately.'
              },
              {
                q: 'Is there a free trial?',
                a: 'All paid plans come with a 14-day free trial. No credit card required to start.'
              },
              {
                q: 'What happens if I exceed my limits?',
                a: 'We\'ll notify you when you\'re approaching limits. You can upgrade anytime or purchase add-ons.'
              },
              {
                q: 'Do you offer academic discounts?',
                a: 'Yes! We offer 50% off for verified academic institutions and students.'
              },
              {
                q: 'Can I get a custom plan?',
                a: 'Enterprise customers can work with us to create custom plans tailored to their needs.'
              },
              {
                q: 'What payment methods do you accept?',
                a: 'We accept all major credit cards, ACH transfers, and wire transfers for enterprise.'
              }
            ].map((faq) => (
              <Grid item xs={12} md={6} key={faq.q}>
                <Paper sx={{ p: 3 }}>
                  <Typography variant="h6" gutterBottom color="primary">
                    {faq.q}
                  </Typography>
                  <Typography variant="body2" color="textSecondary">
                    {faq.a}
                  </Typography>
                </Paper>
              </Grid>
            ))}
          </Grid>
        </Box>

        {/* CTA */}
        <Paper sx={{ 
          p: 6, 
          textAlign: 'center',
          background: 'linear-gradient(135deg, #00e5ff20 0%, #76ff0320 100%)'
        }}>
          <Typography variant="h4" gutterBottom>
            Ready to Accelerate Your Drug Discovery?
          </Typography>
          <Typography variant="body1" color="textSecondary" sx={{ mb: 3 }}>
            Join thousands of researchers using CriOS to discover breakthrough medicines
          </Typography>
          <Button
            variant="contained"
            size="large"
            href="/signup"
            sx={{
              background: 'linear-gradient(45deg, #00e5ff 30%, #76ff03 90%)',
              px: 4,
              py: 1.5
            }}
          >
            Start Your Free Trial
          </Button>
        </Paper>
      </Container>
    </ThemeProvider>
  );
}