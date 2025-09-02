"""
Enterprise Licensing System for CriOS
Handles trial licenses, wire transfers, and enterprise agreements
"""

from datetime import datetime, timedelta
from typing import Dict, Optional, List, Tuple
from enum import Enum
from pydantic import BaseModel, Field
import hashlib
import secrets
import json
from decimal import Decimal

class LicenseType(str, Enum):
    """License tiers with associated pricing"""
    POC = "proof_of_concept"           # $50,000 - 90 days
    PILOT = "pilot_program"             # $99,000 - 6 months
    COMMERCIAL = "commercial"           # $1,500,000 - Perpetual
    ENTERPRISE = "enterprise"           # $5,000,000 - Perpetual
    GLOBAL = "global_transfer"          # $15,000,000 - Full IP transfer

class PaymentStatus(str, Enum):
    """Payment tracking states"""
    PENDING = "pending"
    DEPOSIT_RECEIVED = "deposit_received"
    PARTIAL = "partial_payment"
    COMPLETED = "completed"
    REFUNDED = "refunded"

class LicenseStatus(str, Enum):
    """License activation states"""
    DRAFT = "draft"
    AWAITING_PAYMENT = "awaiting_payment"
    ACTIVE = "active"
    EXPIRED = "expired"
    SUSPENDED = "suspended"
    TERMINATED = "terminated"

# ============================================================================
# LICENSE CONFIGURATION
# ============================================================================

LICENSE_CONFIG = {
    LicenseType.POC: {
        "price": 50000,
        "duration_days": 90,
        "agents": 50,
        "compounds_limit": 100000,
        "pipelines_limit": 3,
        "users_limit": 5,
        "commercial_use": False,
        "support_level": "business_hours",
        "credit_applicable": True,
        "credit_window_days": 120,
    },
    LicenseType.PILOT: {
        "price": 99000,
        "duration_days": 180,
        "agents": 194,
        "compounds_limit": 1000000,
        "pipelines_limit": -1,  # Unlimited
        "users_limit": 20,
        "commercial_use": True,
        "commercial_revenue_cap": 500000,
        "support_level": "24/7",
        "credit_applicable": True,
        "credit_window_days": 270,
        "includes_training": True,
    },
    LicenseType.COMMERCIAL: {
        "price": 1500000,
        "duration_days": -1,  # Perpetual
        "agents": 194,
        "compounds_limit": -1,
        "pipelines_limit": -1,
        "users_limit": -1,
        "commercial_use": True,
        "support_years": 1,
        "source_code_escrow": True,
        "training_users": 50,
    },
    LicenseType.ENTERPRISE: {
        "price": 5000000,
        "duration_days": -1,
        "agents": 194,
        "compounds_limit": -1,
        "pipelines_limit": -1,
        "users_limit": -1,
        "commercial_use": True,
        "white_label": True,
        "custom_agents": 5,
        "support_years": 3,
        "dedicated_support": True,
        "sublicense_rights": True,
    },
    LicenseType.GLOBAL: {
        "price": 15000000,
        "duration_days": -1,
        "full_source_code": True,
        "ip_transfer": True,
        "joint_development": True,
        "support_years": 5,
        "patent_rights": True,
        "unlimited_everything": True,
    }
}

# ============================================================================
# DATA MODELS
# ============================================================================

class LicenseKey(BaseModel):
    """Cryptographically secure license key"""
    key: str
    fingerprint: str
    created_at: datetime
    hardware_id: Optional[str] = None
    
    @classmethod
    def generate(cls, org_id: str, license_type: str) -> "LicenseKey":
        """Generate unique license key"""
        timestamp = datetime.now().isoformat()
        raw_key = f"{org_id}-{license_type}-{timestamp}-{secrets.token_hex(16)}"
        
        # Create key components
        prefix = "CRIOS"
        type_code = license_type[:3].upper()
        unique_id = secrets.token_hex(8).upper()
        checksum = hashlib.sha256(raw_key.encode()).hexdigest()[:8].upper()
        
        # Format: CRIOS-POC-XXXX-XXXX-XXXX-XXXX
        key = f"{prefix}-{type_code}-{unique_id[:4]}-{unique_id[4:8]}-{unique_id[8:12]}-{checksum}"
        
        return cls(
            key=key,
            fingerprint=hashlib.sha256(key.encode()).hexdigest(),
            created_at=datetime.now()
        )

class WireTransferDetails(BaseModel):
    """Wire transfer instructions"""
    bank_name: str = "JPMorgan Chase Bank, N.A."
    bank_address: str = "270 Park Avenue, New York, NY 10017, USA"
    swift_code: str = "CHASUS33"
    routing_number: str = "021000021"
    account_number: str = "XXXXXX1234"  # Masked for security
    account_name: str = "Crowe Research Inc."
    reference: str
    amount: Decimal
    currency: str = "USD"
    
    def to_instructions(self) -> str:
        """Generate formatted wire instructions"""
        return f"""
WIRE TRANSFER INSTRUCTIONS
==========================
Beneficiary: {self.account_name}
Bank: {self.bank_name}
Address: {self.bank_address}
SWIFT/BIC: {self.swift_code}
Routing (ABA): {self.routing_number}
Account: {self.account_number}
Reference: {self.reference}
Amount: {self.currency} {self.amount:,.2f}

IMPORTANT: Include reference number in wire details.
Confirmation email: payments@crios.ai
"""

class LicenseAgreement(BaseModel):
    """Enterprise license agreement"""
    agreement_id: str
    organization: str
    organization_id: str
    license_type: LicenseType
    license_key: Optional[LicenseKey] = None
    
    # Contacts
    primary_contact: str
    primary_email: str
    technical_contact: Optional[str] = None
    billing_contact: Optional[str] = None
    
    # Terms
    start_date: datetime
    end_date: Optional[datetime] = None
    auto_renew: bool = False
    
    # Payment
    total_amount: Decimal
    payment_terms: str = "Net 30"
    payment_status: PaymentStatus = PaymentStatus.PENDING
    payment_milestones: List[Dict] = Field(default_factory=list)
    
    # Status
    status: LicenseStatus = LicenseStatus.DRAFT
    executed_date: Optional[datetime] = None
    signed_by: Optional[str] = None
    
    # Credits from trials
    applied_credit: Decimal = Decimal("0")
    credit_from_license: Optional[str] = None
    
    # Usage tracking
    usage_stats: Dict = Field(default_factory=dict)
    last_activity: Optional[datetime] = None

# ============================================================================
# LICENSE MANAGER
# ============================================================================

class LicenseManager:
    """Manage enterprise licenses and payments"""
    
    def __init__(self):
        self.agreements: Dict[str, LicenseAgreement] = {}
        self.payment_processor = PaymentProcessor()
    
    def create_trial_license(
        self,
        organization: str,
        contact_email: str,
        license_type: LicenseType,
        apply_credit_from: Optional[str] = None
    ) -> Tuple[LicenseAgreement, WireTransferDetails]:
        """Create trial license (POC or Pilot)"""
        
        if license_type not in [LicenseType.POC, LicenseType.PILOT]:
            raise ValueError("Only POC and PILOT licenses can be created as trials")
        
        config = LICENSE_CONFIG[license_type]
        agreement_id = f"AGR-{secrets.token_hex(8).upper()}"
        org_id = hashlib.md5(organization.encode()).hexdigest()[:8].upper()
        
        # Calculate pricing with credit
        total_amount = Decimal(str(config["price"]))
        applied_credit = Decimal("0")
        
        if apply_credit_from:
            # Check if previous license exists and is eligible for credit
            if apply_credit_from in self.agreements:
                prev_agreement = self.agreements[apply_credit_from]
                if prev_agreement.license_type == LicenseType.POC:
                    # Check if within credit window
                    credit_window = timedelta(days=LICENSE_CONFIG[LicenseType.POC]["credit_window_days"])
                    if datetime.now() - prev_agreement.start_date <= credit_window:
                        applied_credit = Decimal(str(LICENSE_CONFIG[LicenseType.POC]["price"]))
                        total_amount -= applied_credit
        
        # Create agreement
        agreement = LicenseAgreement(
            agreement_id=agreement_id,
            organization=organization,
            organization_id=org_id,
            license_type=license_type,
            primary_contact=contact_email.split('@')[0],
            primary_email=contact_email,
            start_date=datetime.now(),
            end_date=datetime.now() + timedelta(days=config["duration_days"]),
            total_amount=total_amount,
            applied_credit=applied_credit,
            credit_from_license=apply_credit_from,
            payment_terms="Net 30",
            status=LicenseStatus.AWAITING_PAYMENT
        )
        
        # Generate license key
        agreement.license_key = LicenseKey.generate(org_id, license_type.value)
        
        # Create wire transfer details
        wire_details = WireTransferDetails(
            reference=agreement_id,
            amount=total_amount,
            currency="USD"
        )
        
        # Store agreement
        self.agreements[agreement_id] = agreement
        
        return agreement, wire_details
    
    def create_enterprise_license(
        self,
        organization: str,
        contact_email: str,
        license_type: LicenseType,
        payment_milestones: Optional[List[Dict]] = None
    ) -> Tuple[LicenseAgreement, List[WireTransferDetails]]:
        """Create enterprise license with milestone payments"""
        
        if license_type not in [LicenseType.COMMERCIAL, LicenseType.ENTERPRISE, LicenseType.GLOBAL]:
            raise ValueError("Not an enterprise license type")
        
        config = LICENSE_CONFIG[license_type]
        agreement_id = f"ENT-{secrets.token_hex(8).upper()}"
        org_id = hashlib.md5(organization.encode()).hexdigest()[:8].upper()
        
        total_amount = Decimal(str(config["price"]))
        
        # Default payment milestones
        if not payment_milestones:
            payment_milestones = [
                {"percentage": 10, "trigger": "signing", "description": "Upon contract execution"},
                {"percentage": 40, "trigger": "delivery", "description": "Upon system delivery"},
                {"percentage": 50, "trigger": "acceptance", "description": "Upon acceptance testing"},
            ]
        
        # Create agreement
        agreement = LicenseAgreement(
            agreement_id=agreement_id,
            organization=organization,
            organization_id=org_id,
            license_type=license_type,
            primary_contact=contact_email.split('@')[0],
            primary_email=contact_email,
            start_date=datetime.now(),
            end_date=None,  # Perpetual
            total_amount=total_amount,
            payment_milestones=payment_milestones,
            payment_terms="Milestone-based",
            status=LicenseStatus.DRAFT
        )
        
        # Generate license key
        agreement.license_key = LicenseKey.generate(org_id, license_type.value)
        
        # Create wire instructions for each milestone
        wire_instructions = []
        for i, milestone in enumerate(payment_milestones):
            milestone_amount = total_amount * Decimal(str(milestone["percentage"] / 100))
            wire = WireTransferDetails(
                reference=f"{agreement_id}-M{i+1}",
                amount=milestone_amount,
                currency="USD"
            )
            wire_instructions.append(wire)
        
        # Store agreement
        self.agreements[agreement_id] = agreement
        
        return agreement, wire_instructions
    
    def activate_license(self, agreement_id: str, payment_confirmation: str) -> bool:
        """Activate license after payment confirmation"""
        
        if agreement_id not in self.agreements:
            return False
        
        agreement = self.agreements[agreement_id]
        
        # Verify payment
        if self.payment_processor.verify_payment(agreement_id, payment_confirmation):
            agreement.status = LicenseStatus.ACTIVE
            agreement.executed_date = datetime.now()
            agreement.payment_status = PaymentStatus.COMPLETED
            
            # Send activation email
            self._send_activation_email(agreement)
            
            return True
        
        return False
    
    def check_license_validity(self, license_key: str) -> Dict:
        """Validate license key and return entitlements"""
        
        for agreement in self.agreements.values():
            if agreement.license_key and agreement.license_key.key == license_key:
                config = LICENSE_CONFIG[agreement.license_type]
                
                # Check expiration for trial licenses
                if agreement.end_date and datetime.now() > agreement.end_date:
                    return {"valid": False, "reason": "License expired"}
                
                if agreement.status != LicenseStatus.ACTIVE:
                    return {"valid": False, "reason": f"License status: {agreement.status}"}
                
                return {
                    "valid": True,
                    "organization": agreement.organization,
                    "type": agreement.license_type.value,
                    "entitlements": {
                        "agents": config.get("agents", 194),
                        "compounds_limit": config.get("compounds_limit", -1),
                        "pipelines_limit": config.get("pipelines_limit", -1),
                        "users_limit": config.get("users_limit", -1),
                        "commercial_use": config.get("commercial_use", False),
                        "expires": agreement.end_date.isoformat() if agreement.end_date else None
                    }
                }
        
        return {"valid": False, "reason": "Invalid license key"}
    
    def _send_activation_email(self, agreement: LicenseAgreement):
        """Send license activation email"""
        # Implementation would use SendGrid or similar
        pass

# ============================================================================
# PAYMENT PROCESSOR
# ============================================================================

class PaymentProcessor:
    """Handle wire transfer verification and tracking"""
    
    def __init__(self):
        self.pending_payments: Dict[str, Dict] = {}
        self.confirmed_payments: Dict[str, Dict] = {}
    
    def verify_payment(self, agreement_id: str, confirmation_code: str) -> bool:
        """Verify wire transfer confirmation"""
        
        # In production, this would:
        # 1. Check with banking API for wire confirmation
        # 2. Verify amount matches agreement
        # 3. Log transaction details
        # 4. Trigger compliance checks for large amounts
        
        # Simplified verification
        if confirmation_code.startswith("WIRE-") and len(confirmation_code) == 20:
            self.confirmed_payments[agreement_id] = {
                "confirmation": confirmation_code,
                "timestamp": datetime.now(),
                "verified": True
            }
            return True
        
        return False
    
    def initiate_refund(self, agreement_id: str, amount: Decimal, reason: str) -> str:
        """Process refund for canceled agreement"""
        
        refund_id = f"REF-{secrets.token_hex(8).upper()}"
        
        # In production, this would:
        # 1. Initiate wire reversal with bank
        # 2. Create audit trail
        # 3. Update accounting systems
        # 4. Notify compliance team
        
        return refund_id

# ============================================================================
# API ENDPOINTS FOR LICENSING
# ============================================================================

def create_licensing_endpoints(app):
    """Add licensing endpoints to FastAPI app"""
    
    license_manager = LicenseManager()
    
    @app.post("/api/license/trial")
    async def create_trial(
        organization: str,
        contact_email: str,
        license_type: str,
        apply_credit_from: Optional[str] = None
    ):
        """Create trial license (POC or Pilot)"""
        
        try:
            license_type_enum = LicenseType(license_type)
            agreement, wire_details = license_manager.create_trial_license(
                organization,
                contact_email,
                license_type_enum,
                apply_credit_from
            )
            
            return {
                "agreement_id": agreement.agreement_id,
                "license_key": agreement.license_key.key,
                "amount": str(agreement.total_amount),
                "applied_credit": str(agreement.applied_credit),
                "wire_instructions": wire_details.to_instructions(),
                "expires": agreement.end_date.isoformat()
            }
        except Exception as e:
            return {"error": str(e)}
    
    @app.post("/api/license/enterprise")
    async def create_enterprise(
        organization: str,
        contact_email: str,
        license_type: str
    ):
        """Create enterprise license"""
        
        try:
            license_type_enum = LicenseType(license_type)
            agreement, wire_instructions = license_manager.create_enterprise_license(
                organization,
                contact_email,
                license_type_enum
            )
            
            return {
                "agreement_id": agreement.agreement_id,
                "license_key": agreement.license_key.key,
                "total_amount": str(agreement.total_amount),
                "payment_milestones": agreement.payment_milestones,
                "wire_instructions": [w.to_instructions() for w in wire_instructions]
            }
        except Exception as e:
            return {"error": str(e)}
    
    @app.post("/api/license/activate")
    async def activate_license(
        agreement_id: str,
        payment_confirmation: str
    ):
        """Activate license after payment"""
        
        if license_manager.activate_license(agreement_id, payment_confirmation):
            return {"status": "activated", "message": "License successfully activated"}
        else:
            return {"status": "failed", "message": "Payment verification failed"}
    
    @app.get("/api/license/validate/{license_key}")
    async def validate_license(license_key: str):
        """Validate license key"""
        return license_manager.check_license_validity(license_key)

# Usage example:
"""
# Create POC license
agreement, wire = license_manager.create_trial_license(
    "Pfizer Inc.",
    "research@pfizer.com",
    LicenseType.POC
)
print(wire.to_instructions())

# Upgrade to Pilot with credit
agreement2, wire2 = license_manager.create_trial_license(
    "Pfizer Inc.",
    "research@pfizer.com", 
    LicenseType.PILOT,
    apply_credit_from=agreement.agreement_id  # Apply $50K credit
)
# Only pay $49K instead of $99K

# Create enterprise license with milestones
enterprise, wires = license_manager.create_enterprise_license(
    "Johnson & Johnson",
    "licensing@jnj.com",
    LicenseType.GLOBAL
)
# Returns 3 wire instructions for $1.5M, $6M, $7.5M milestones
"""